package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Bytes;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple3;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary        = "Tool to genotype inversions on a single sample, based on assembled contigs from previous steps in the SV pipeline, " +
                         "and outputs genotype likelihoods.",
        oneLineSummary = "Genotype single sample inversion variant on Spark.",
        programGroup   = StructuralVariationSparkProgramGroup.class)
public final class GenotypingSingleSampleInversionSpark extends GenotypingSingleSampleSVSpark {
    private static final long serialVersionUID = 1L;

    private static final int readLength = 151;
    private static final double expectedBaseErrorRate = 1.0/readLength; // blunt guess, will change
    // LikelihoodEngineArgumentCollection.phredScaledGlobalReadMismappingRate default value, used in HC
    private static final double maximumLikelihoodDifferenceCap = 45;

    // -----------------------------------------------------------------------------------------------
    // Overrides
    // -----------------------------------------------------------------------------------------------

    @Override
    public void runTool(final JavaSparkContext ctx){

        final InversionReadLikelihoodCalculator rllCalculator = new InversionReadLikelihoodCalculator();
        rllCalculator.initialize();
        rllCalculator.configure();

        super.readLikelihoodCalculator = rllCalculator;
        super.genotypeLikelihoodCalculator = new GenotypeLikelihoodCalculators().getInstance(2, 2); // TODO: diploid, biallelic assumption

        super.runTool(ctx);
    }

    /**
     * Custom class for holding information around two inversion breakpoints.
     */
    @VisibleForTesting
    static final class InversionJunction extends SVJunction {
        private static final long serialVersionUID = 1L;

        private final BreakpointAllele.InversionType invType;


        InversionJunction(final VariantContext vc,
                          final ReferenceMultiSource reference,
                          final Map<Long, List<LocalAssemblyContig>> assemblyID2AlignmentRecords){

            super(vc, assemblyID2AlignmentRecords);

            final int isFiveToThree = vc.getAttributeAsBoolean(GATKSVVCFHeaderLines.INV_5_TO_3, false) ? 0b1 : 0b0;
            final int bit = isFiveToThree + (vc.getAttributeAsBoolean(GATKSVVCFHeaderLines.INV_3_TO_5, false) ? 0b10 : 0b0);
            if (bit == 0b0) {
                invType = BreakpointAllele.InversionType.INV_NONE;
            } else if (bit == 0b01) {
                invType = BreakpointAllele.InversionType.INV_5_TO_3;
            } else if (bit == 0b10) {
                invType = BreakpointAllele.InversionType.INV_3_TO_5;
            } else {
                throw new IllegalArgumentException("Seemingly broken VCF, where the inversion breakpoint is of both type 5-to-3 and 3-to-5. Site: "
                        + vc.getContig() + ":" + vc.getID());
            }
            setWhichEnd();
            constructAlleles(reference);
        }

        @VisibleForTesting
        @Override
        protected void setWhichEnd(){
            if (invType == BreakpointAllele.InversionType.INV_NONE) {
                whichEnd = 0;
            } else if (invType == BreakpointAllele.InversionType.INV_5_TO_3) {
                whichEnd = 1;
            } else if (invType == BreakpointAllele.InversionType.INV_3_TO_5) {
                whichEnd = 2;
            } else {
                throw new IllegalArgumentException("Seemingly broken VCF, un-recognizable inversion type, Site: "
                        + vc.getContig() + vc.getStart());
            }
        }

        /**
         * TODO: get the short window case correct.
         */
        @VisibleForTesting
        @Override
        protected Tuple3<String, String, int[]> constructReferenceWindows(){

            final List<SimpleInterval> leftAlignedBreakpointLocations = getLeftAlignedBreakpointLocations();
            final SimpleInterval fiveEndBPLoc = leftAlignedBreakpointLocations.get(0);
            final SimpleInterval threeEndBPLoc = leftAlignedBreakpointLocations.get(1);

            final String contig = fiveEndBPLoc.getContig();

            final int flank = readLength - 1; // -1, think about it
            int ll, lr, rl, rr;
            if (invType!=BreakpointAllele.InversionType.INV_NONE) {
                final int extraFlanking = (maybeNullHomology==null) ? 0 : maybeNullHomology.length;
                ll = fiveEndBPLoc.getStart()  - flank;
                lr = fiveEndBPLoc.getEnd()    + flank + extraFlanking;
                rl = threeEndBPLoc.getStart() - flank - extraFlanking;
                rr = threeEndBPLoc.getEnd()   + flank;
            } else {
                ll = fiveEndBPLoc.getStart()  - flank;
                lr = fiveEndBPLoc.getEnd()    + flank;
                rl = threeEndBPLoc.getStart() - flank;
                rr = threeEndBPLoc.getEnd()   + flank;
            }
            return new Tuple3<>(contig, contig, new int[]{ll, lr, rl, rr});
        }

        @Override
        protected final List<SVDummyAllele> constructAlternateAlleles(final ReferenceMultiSource reference){
            // TODO: because the assembler may decide to stop short of extending to the full windows around identified breakpoint,
            //       making use of the assembled contig to construct the alt allele is not mature yet until we can make sense out of the assembly graph
            //       so now we simply use the reference bases to construct the alt allele,
            //       what could be done intermediately is to amend the allele bases with the contigs, although this way
            //       there will be more one one alt allele for one end, because that's the reason why the assembler decides keep the bubble there

            try{
                byte[] tempBases1 = reference.getReferenceBases(null, new SimpleInterval(referenceWindows._1(), referenceWindows._3()[0], fiveEndBPLoc.getStart())).getBases();
                byte[] tempBases2 = reference.getReferenceBases(null, new SimpleInterval(referenceWindows._1(), referenceWindows._3()[2], threeEndBPLoc.getStart())).getBases();
                if(maybeNullInsertedSeq!=null){
                    tempBases1 = Bytes.concat(tempBases1, maybeNullInsertedSeq);
                }
                SequenceUtil.reverseComplement(tempBases2);
                final SVDummyAllele fiveToThreeAltAllele = new SVDummyAllele(Bytes.concat(tempBases1, tempBases2), false);

                tempBases1 = reference.getReferenceBases(null, new SimpleInterval(referenceWindows._1(), fiveEndBPLoc.getStart(), referenceWindows._3()[1])).getBases();
                tempBases2 = reference.getReferenceBases(null, new SimpleInterval(referenceWindows._1(), threeEndBPLoc.getStart(), referenceWindows._3()[3])).getBases();
                SequenceUtil.reverseComplement(tempBases1);
                if(maybeNullInsertedSeq!=null){
                    tempBases2 = Bytes.concat(maybeNullInsertedSeq, tempBases2);
                }
                final SVDummyAllele threeToFiveAltAllele = new SVDummyAllele(Bytes.concat(tempBases1, tempBases2), false);

                return Arrays.asList(fiveToThreeAltAllele, threeToFiveAltAllele);
            } catch (final IOException ioex){
                throw new GATKException("Cannot resolve reference for constructing ref allele for inversion junction");
            }
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            if (!super.equals(o)) return false;

            final InversionJunction that = (InversionJunction) o;

            return invType == that.invType;

        }

        @Override
        public int hashCode() {
            int result = super.hashCode();
            result = 31 * result + invType.ordinal();
            return result;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected InversionJunction convertToSVJunction(final VariantContext vc,
                                                    final Broadcast<Map<Long, List<LocalAssemblyContig>>> assembly2Alignments,
                                                    final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast){
        return new InversionJunction(vc, referenceMultiSourceBroadcast.getValue(), assembly2Alignments.getValue());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean readSuitableForGenotypingJunction(final SVJunction junction, final GATKRead read){
        return readResideInJunctionWindow((InversionJunction) junction, read); // TODO: read could also share k-mer with ref or alt
    }

    /**
     * {@inheritDoc}
     *
     * Normalizes the read likelihoods, see {@link ReadLikelihoods#normalizeLikelihoods(boolean, double)}.
     * filter out poorly modeled reads, see {@link ReadLikelihoods#filterPoorlyModeledReads(double)},
     * and
     * marginalizes the read likelihoods by reducing from multiple ref alleles and multiple alt alleles
     * used in the RLC (read likelihood calculation) step to a single symbolic ref allele and
     * a single symbolic alt allele (see {@link ReadLikelihoods#marginalize(Map)} for logic)
     * so that the following assumption is made:
     * (P=2, A=2) =&gt; [0/0, 0/1, 1/1] three possible genotypes for a diploid biallelic sample.
     * TODO: is this marginalization logic correct?
     *
     * TODO: Re-alignment of reads back to their best allele is not implemented yet (as it happens in HC).
     */
    @Override
    protected ReadLikelihoods<SVDummyAllele> updateReads(final ReadLikelihoods<SVDummyAllele> matrix,
                                                         final SVJunction inversionJunction){

        matrix.normalizeLikelihoods(false, QualityUtils.qualToErrorProbLog10(maximumLikelihoodDifferenceCap));
        matrix.filterPoorlyModeledReads(expectedBaseErrorRate);

        final List<SVDummyAllele> alleles = inversionJunction.getAlleles();
        final List<SVDummyAllele> vcAlleles = inversionJunction.getOriginalVC().getAlleles().stream().map(SVDummyAllele::new).collect(Collectors.toList());

        final Map<SVDummyAllele, List<SVDummyAllele>> symbolicOnes2RealAlleles = new HashMap<>();
        symbolicOnes2RealAlleles.put(vcAlleles.get(0), alleles.stream().filter(SVDummyAllele::isReference).collect(Collectors.toList()));
        symbolicOnes2RealAlleles.put(vcAlleles.get(1), alleles.stream().filter(SVDummyAllele::isNonReference).collect(Collectors.toList()));
        return matrix.marginalize(symbolicOnes2RealAlleles);
    }

    // -----------------------------------------------------------------------------------------------
    // Utilities
    // -----------------------------------------------------------------------------------------------

    /**
     * Test if a particular read starts after or ends before the two windows spanned by the inversion junction.
     */
    @VisibleForTesting
    boolean readResideInJunctionWindow(final InversionJunction invJunction, final GATKRead read){

        if(read.isUnmapped()) return false; // TODO: really stupid!

        final Tuple3<String, String, int[]> windows = invJunction.getReferenceWindows();
        final String contig = windows._1();
        final SimpleInterval left = new SimpleInterval(contig, windows._3()[0], windows._3()[1]);
        final SimpleInterval right = new SimpleInterval(contig, windows._3()[2], windows._3()[3]);

        if ( !contig.equalsIgnoreCase(read.getContig())
                || ! (left.overlaps(read) || right.overlaps(read)) ) {
            return false;
        } else { // read start after certain locations or end before certain locations
            final int[] boundaries = windows._3();
            final int readUnclippedStart = read.getUnclippedStart();
            final int readUnclippedEnd = read.getUnclippedEnd();

            return readUnclippedStart>=boundaries[0] || readUnclippedStart>=boundaries[2] ||
                   readUnclippedEnd<=boundaries[1] || readUnclippedEnd<=boundaries[3];
        }
    }
}
