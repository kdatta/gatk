package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import scala.Tuple2;
import scala.Tuple3;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * TODO: when we move to genotyping multiple types of sv, this should be moved to an engine class instead of base tool class.
 * TODO: two major methods: getBamReads() and getFASTQReads() are difficult to test.
 * TODO: right now we are dealing with single sample only, so we don't need a {@link org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingData} alike,
 *       for the purpose of storing a {@link org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel}
 *       and getting the appropriate {@link org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix} for the sample.
 *       it may make sense to refactor {@link ReadLikelihoods} for a more explicit sample-stratified structure
 */
abstract class GenotypingSingleSampleSVSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            optional  = false)
    protected String VCFFile;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = "fastq",
            fullName  = "FASTQDir",
            optional  = false)
    protected String pathToFASTQFiles;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = "assembly",
            fullName  = "assembledFASTADir",
            optional  = false)
    protected String pathToAssembledFASTAFiles;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = "aln",
            fullName  = "contigAlignments",
            optional  = false)
    protected String pathToContigAlignments;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    protected String outFileName;


    @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
            fullName = "fastaReference", optional = false)
    private String fastaReference;

    @Argument(doc = "Path to a file where debugging string will be writen out to", shortName = "dso",
            fullName = "debugstringoutpath", optional = true)
    @Advanced
    private String debugStringSavingOutputPath;

    // read likelihood calculator and subsequent single sample genotyper
    protected SVReadLikelihoodCalculator readLikelihoodCalculator;
    protected GenotypeLikelihoodCalculator genotypeLikelihoodCalculator;

    @Override
    public final boolean requiresReference() {
        return true;
    }

    @Override
    public final boolean requiresReads() {
        return true;
    }

    // for developer performance debugging use
    static final boolean in_debug_state = true;
    private static final Logger logger = LogManager.getLogger(GenotypingSingleSampleSVSpark.class);

    static final String testSampleName = "NA12878_PCR-_30X";

    // -----------------------------------------------------------------------------------------------
    // Manager
    // -----------------------------------------------------------------------------------------------

    @Override
    public void runTool(final JavaSparkContext ctx){

        // parse discovery VC calls, gather and filter reads that will be used in genotyping
        final List<SVJunction> junctions = parseVCFInput(ctx, VCFFile);
        filterJunctions(junctions, makeJunctionFilter());

        final JavaPairRDD<SVJunction, List<GATKRead>> genotypingReads = prepareReads(ctx, junctions, pathToFASTQFiles).cache();
        final Set<SVJunction> sparkJunctions = genotypingReads.keys().collect().stream().collect(Collectors.toSet());
        if(Sets.symmetricDifference(sparkJunctions, junctions.stream().collect(Collectors.toSet())).size()!=0){
            throw new IllegalStateException("Junctions are different before and after read gathering!");
        }
        logger.debug(".........DONE MERGING READS.........");

        // send reads to read likelihood calculator, but first configure the calculator appropriately
        // TODO: think about how to configure calculator
        readLikelihoodCalculator.initialize();
        readLikelihoodCalculator.configure();

        // pre-process data like HC (or rather PairHMMLikelihoodCalculationEngine) does
        final JavaPairRDD<SVJunction, List<GATKRead>> preprocessedGenotypingReads = genotypingReads.mapValues(reads -> readLikelihoodCalculator.preprocessReads(reads)).cache();
        genotypingReads.unpersist();
        { preprocessedGenotypingReads.count();logger.debug(".........DONE PREPROCESSING READS........."); }

        // compute read likelihoods
        final JavaPairRDD<SVJunction, ReadLikelihoods<SVDummyAllele>> step1 = preprocessedGenotypingReads.mapToPair(this::getSvJunctionReadLikelihoods).cache();
//        if(in_debug_state) writeoutDebugReadLikelihoods(step1, debugStringSavingOutputPath);
        { step1.count();logger.debug(".........RLC DONE........."); }
        preprocessedGenotypingReads.unpersist();
        readLikelihoodCalculator.close();

        // post processing of read likelihoods, including normalization and (optional) re-alignment
        final JavaPairRDD<SVJunction, ReadLikelihoods<SVDummyAllele>> step2 = step1.mapToPair(p -> new Tuple2<>(p._1(), updateReads(p._2(), p._1()))).cache();
        { step2.count();logger.debug(".........POSTPROCESSING DONE........."); }
        step1.unpersist();

        // compute single sample genotype likelihoods
        final JavaPairRDD<SVJunction, GenotypeLikelihoods> step3 = step2.mapValues(gl -> genotypeLikelihoodCalculator.genotypeLikelihoods(  gl.sampleMatrix(0)  )).cache();
        { step3.count();logger.debug(".........GLC DONE........."); }
        step2.unpersist();

        // turn back to VC and save
        final SAMSequenceDictionary referenceSequenceDictionary = new ReferenceMultiSource(getAuthenticatedGCSOptions(), fastaReference, ReferenceWindowFunctions.IDENTITY_FUNCTION).getReferenceSequenceDictionary(null);
        final List<VariantContext> genotypedVCs = collectAndSortGenotypedVC(step3, referenceSequenceDictionary);
        output(genotypedVCs, referenceSequenceDictionary, getAuthenticatedGCSOptions(), outFileName);
    }

    // -----------------------------------------------------------------------------------------------
    // To be overridden
    // -----------------------------------------------------------------------------------------------

    /**
     * Converting from a discovery {@link VariantContext} call to an appropriate {@link SVJunction} suitable for genotyping.
     */
    protected abstract SVJunction convertToSVJunction(final VariantContext vc,
                                                      final Broadcast<Map<Long, List<LocalAssemblyContig>>> assemblyID2assemblyContents,
                                                      final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast);

    @VisibleForTesting
    protected Predicate<SVJunction> makeJunctionFilter(){
        return new JunctionFilterBasedOnInterval();
    }

    /**
     * Go back to the original bam and gather reads around region identified by {@link SVJunction}.
     *
     * Reads for each junction are filtered based on subclass custom defined {@link #readSuitableForGenotypingJunction(SVJunction, GATKRead)}.
     *
     * Note that read names are changed by appending "/1" or "/2" to their original names so later when bam reads are
     * merged with FASTQ reads extracted by {@link #getFASTQReads(JavaSparkContext, String, List)},
     * duplicate reads can be removed by read names.
     */
    protected JavaPairRDD<SVJunction, List<GATKRead>> getBamReads(final List<SVJunction> svJunctions){
        return getReads()
                .filter(read -> svJunctions.stream().anyMatch(junction -> readSuitableForGenotypingJunction(junction, read))) // pre-filer: if suitable for any junction
                .flatMapToPair( read -> svJunctions.stream()
                                                    .filter(junction -> readSuitableForGenotypingJunction(junction, read))
                                                    .map(junc -> new Tuple2<>(junc, read))
                                                    .collect(Collectors.toList()) )                                           // read -> { (junction, read) }
                .groupByKey()
                .mapValues( it -> {// bwa strips out the "/1" "/2" part, but we need these for merging with FASTQ reads
                            it.forEach(read -> read.setName(read.getName() + (read.isFirstOfPair() ? "/1" : "/2")));
                            return StreamSupport.stream(it.spliterator(), false).collect(Collectors.toList());}); // turn iterable to list
    }

    /**
     * Filter read that will be used for genotyping.
     * Assumes the logic applies to all junctions of the same type.
     * TODO: this interface implicitly works against pairedness.
     */
    protected abstract boolean readSuitableForGenotypingJunction(final SVJunction junction,
                                                                 final GATKRead read);

    /**
     * Update read based on result from its realignment to ref and alt contigs.
     * Input reads are deep copied for making the necessary changes.
     */
    protected abstract ReadLikelihoods<SVDummyAllele> updateReads(final ReadLikelihoods<SVDummyAllele> reads,
                                                                  final SVJunction junction);

    // -----------------------------------------------------------------------------------------------
    // Utilities: delegation (to avoid lambda calls that requires too many other classes to be serialized)
    // -----------------------------------------------------------------------------------------------

    // TODO: dummy place holder for a list of samples, subject to future change, which requires some of the managed structs and interface and implementation to change
    private Tuple2<SVJunction, ReadLikelihoods<SVDummyAllele>> getSvJunctionReadLikelihoods(final Tuple2<SVJunction, List<GATKRead>> p) {
        final SVJunction junction = p._1();
        final Map<String, List<GATKRead>> sample2Reads = new LinkedHashMap<>();

        final SampleList sampleList = SampleList.singletonSampleList(testSampleName);
        sample2Reads.put(sampleList.getSample(0),p._2());

        return new Tuple2<>(junction, readLikelihoodCalculator.computeReadLikelihoods(sampleList, p._2(), sample2Reads, junction));
    }

    /**
     * TODO: this is possible only because of single sample; full-blown cohort genotyping for SV is to be done (SNP and indel models should be learned)
     * Simply infer genotype of the single diploid sample from the computed PL.
     * @return a new genotype based on the PL of input genotype
     */
    @VisibleForTesting
    static Genotype inferGenotypeFromPL(final Genotype gtWithPL, final List<Allele> alleles){

        final GenotypeBuilder builder = new GenotypeBuilder(gtWithPL);
        if(!GATKVariantContextUtils.isInformative(gtWithPL.getLikelihoods().getAsVector())){
            builder.noPL();
        }

        final double[] ll = MathUtils.normalizeFromLog10(gtWithPL.getLikelihoods().getAsVector(), false, true);

        GATKVariantContextUtils.makeGenotypeCall(2, builder, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, ll, alleles);

        return builder.make();
    }

    // -----------------------------------------------------------------------------------------------
    // Utilities: IO
    // -----------------------------------------------------------------------------------------------

    /**
     * Parse input (output from previous step in the pipeline),
     * and return a list of convenience struct for use by specific genotyping module.
     *
     * Depends on subclasses' {@link #convertToSVJunction(VariantContext, Broadcast, Broadcast)} for the actual mapping.
     */
    private List<SVJunction> parseVCFInput(final JavaSparkContext ctx,
                                           final String vcfInput){

        final List<VariantContext> listOfVCs = new VariantsSparkSource(ctx).getParallelVariantContexts(vcfInput, null).collect();

        final Map<Long, List<LocalAssemblyContig>> assemblyID2assembleContents = mapAssemblyID2ItsAlignments(ctx, pathToAssembledFASTAFiles, pathToContigAlignments, gatherAsmIDs(listOfVCs));
        final Broadcast<Map<Long, List<LocalAssemblyContig>>> assemblyID2assembleContentsBroadcast = ctx.broadcast(assemblyID2assembleContents);

        final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast = ctx.broadcast(getReference());

        if(in_debug_state) logger.info(".........BEGIN CONSTRUCTING JUNCTIONS FROM VC AND ALIGNMENTS.........");
        final List<SVJunction> junctions = listOfVCs.stream().map(vc -> convertToSVJunction(vc, assemblyID2assembleContentsBroadcast, referenceMultiSourceBroadcast)).collect(Collectors.toList());
        if(in_debug_state) logger.info(String.format(".........CONSTRUCTED %d JUNCTIONS FROM VC AND ALIGNMENTS.........", junctions.size()));
        assemblyID2assembleContentsBroadcast.destroy();referenceMultiSourceBroadcast.destroy();

        return junctions;
    }

    @VisibleForTesting
    final class JunctionFilterBasedOnInterval implements Predicate<SVJunction>{
        private List<SimpleInterval> regionsToExclude = Arrays.asList(new SimpleInterval("2", 33141200, 33141700),
                new SimpleInterval("MT", 12867, 14977),
                new SimpleInterval("17", 41379304, 41401057),
                new SimpleInterval("17", 41381590, 41463575),
                new SimpleInterval("17", 41401193, 41466836));

        public boolean test(final SVJunction j){
            final Tuple3<String, String, int[]> windows = j.getReferenceWindows();
            final String contig = windows._1();
            final SimpleInterval left = new SimpleInterval(contig, windows._3()[0], windows._3()[1]);
            final SimpleInterval right = new SimpleInterval(contig, windows._3()[2], windows._3()[3]);
            return regionsToExclude.stream().anyMatch(region -> region.overlaps(left) || region.overlaps(right));
        }
    }

    @VisibleForTesting
    static void filterJunctions(final List<SVJunction> junctions, Predicate<SVJunction> predicate){
        final List<SVJunction> junctionsToBeFiltered = junctions.stream().filter(predicate::test).collect(Collectors.toList());
        if(junctionsToBeFiltered.size()!=0){
            logger.warn(".........FOLLOWING JUNCTIONS TO BE EXCLUDED IN GENOTYPING BECAUSE OF HIGH COVERAGE.........");
            junctionsToBeFiltered.stream().map(j -> j.getOriginalVC().getID()).forEach(logger::warn);
            logger.warn(String.format(".........%d JUNCTIONS LEFT TO BE GENOTYPED.........", junctions.size()));
            junctions.removeAll(junctionsToBeFiltered);
        }
    }

    /**
     * @return a mapping from assembly id to the result of locally assembled contigs aligned back to ref,
     *         given the path to re-alignment result.
     */
    @VisibleForTesting
    static Map<Long, List<LocalAssemblyContig>> mapAssemblyID2ItsAlignments(final JavaSparkContext ctx,
                                                                            final String localAssemblyResultPath,
                                                                            final String alignmentResultPath,
                                                                            final Set<Long> assembliesToKeep){

        final Map<Long, ContigsCollection> asmId2Contigs = ContigsCollection.loadContigsCollectionKeyedByAssemblyId(ctx, localAssemblyResultPath)
                .mapToPair( p -> new Tuple2<>(Long.parseLong(p._1().replace(SVConstants.FASTQ_OUT_PREFIX, "")), p._2()))
                .collectAsMap();

        final JavaPairRDD<Long, Iterable<AlignmentRegion>> asmId2Alignments = ctx.textFile(alignmentResultPath)
                .map(ContigsCollection::parseAlignedAssembledContigLine)
                .mapToPair(ar -> new Tuple2<>(Long.parseLong(ar.assemblyId.replace(SVConstants.FASTQ_OUT_PREFIX, "")), ar))
                .groupByKey();

        // asmId, mapped to a list of LocalAssemblyContig's
        return asmId2Alignments.mapToPair(pair -> {
            final Long asmId = pair._1();
            final Map<String, String> contigId2Seq = asmId2Contigs.get(asmId).getContents().stream().collect(Collectors.toMap(contig -> contig.contigID, contig -> contig.seq));

            final Map<String, List<AlignmentRegion>> contigId2Alignments = StreamSupport.stream(pair._2().spliterator(), false).collect(Collectors.groupingBy(ar -> ar.contigId));

            return new Tuple2<>(asmId, contigId2Seq.entrySet().stream().map(p -> new LocalAssemblyContig(asmId, p.getKey(), p.getValue(), contigId2Alignments.get(p.getKey()))).collect(Collectors.toList()));
        }).filter(p->assembliesToKeep.contains(p._1())).collectAsMap();
    }

    @VisibleForTesting
    static Set<Long> gatherAsmIDs(final List<VariantContext> vcs){
        return vcs.stream()
                .flatMap(vc -> Arrays.stream(vc.getAttributeAsString(GATKSVVCFHeaderLines.ASSEMBLY_IDS, "").replace("[", "").replace("]", "").split(", ")).map(Long::valueOf)  )
                .collect(Collectors.toSet());
    }

    /**
     * Prepare reads that will be used for likelihood calculations.
     *
     * Reads are gathered from two resources:
     * 1) original FASTQ records pulled out by {@link FindBreakpointEvidenceSpark} that are associated with
     *      putative breakpoints, which in turn the appropriate {@link SVJunction} are associated with.
     * 2) original bam file. Exactly how reads are pulled in from the bam file is determined by {@link #getBamReads(List)},
     *      whose exact behavior is overridable by subclasses.
     *
     * For a junction, reads that are gathered from both sources (bam, and reconstructed from FASTQ), will be de-duplicated
     * such that the bam version will be kept.
     * @param svJunctions           a list of sv junctions to be genotyped
     * @param pathToFASTQFiles      a string to the directory where the FASTQ files generated by {@link FindBreakpointEvidenceSpark} are stored
     */
    private JavaPairRDD<SVJunction, List<GATKRead>> prepareReads(final JavaSparkContext ctx,
                                                                 final List<SVJunction> svJunctions,
                                                                 final String pathToFASTQFiles){

        // gather reads from FASTQ files collected in previous steps in pipeline
        final JavaPairRDD<SVJunction, List<GATKRead>> fastqReads = getFASTQReads(ctx, pathToFASTQFiles, svJunctions);

        // gather reads from bam file (will overlap with FASTQ reads)
        if(in_debug_state) logger.info(".........BEGIN EXTRACTING BAM READS.........");
        final JavaPairRDD<SVJunction, List<GATKRead>> bamReads = getBamReads(svJunctions).cache();
        if(in_debug_state){
            bamReads.count(); logger.info(".........DONE EXTRACTING BAM READS.........");
        }

        // merge together the fastqReads and bamReads (deduplication to follow)
        if(in_debug_state) logger.info(".........BEGIN MERGING READS.........");
        final JavaPairRDD<SVJunction, Tuple2<List<GATKRead>, List<GATKRead>>> redundantReads = fastqReads.join(bamReads).cache();
        bamReads.unpersist();

        return redundantReads.mapValues(pair -> mergeReadsByName(pair._1, pair._2));
    }

    /**
     * Get FASTQ reads that were pulled out by {@link FindBreakpointEvidenceSpark}.
     * Reads are filtered by custom filters defined by sub classes (see {@link #readSuitableForGenotypingJunction(SVJunction, GATKRead) <ReferenceMultiSource>)}.
     * @return Reads associated with a particular breakpoint.
     */
    private JavaPairRDD<SVJunction, List<GATKRead>> getFASTQReads(final JavaSparkContext ctx,
                                                                  final String pathToFASTQFiles,
                                                                  final List<SVJunction> svJunctions){

        final Map<Long, List<SVJunction>> asmid2SVJunctions = getAssemblyToJunctionsMapping(svJunctions);
        final Set<Long> interestingAsmIDs = new HashSet<>( asmid2SVJunctions.keySet() );
        logger.debug(String.format(".........%d JUNCTIONS, %d ASSEMBLIES, WAITING TO GATHER READS.........",
                                    StreamSupport.stream(asmid2SVJunctions.values().spliterator(), false).flatMap(List::stream).distinct().collect(Collectors.toList()).size(),
                                    interestingAsmIDs.size()));

        // map assembly id (only those giving signal) to FASTQ contents then reconstruct read
        final JavaPairRDD<Long, List<GATKRead>> preFilteredReads = SVFastqUtils.loadFASTQFiles(ctx, pathToFASTQFiles)
                                                                                .mapToPair(pair -> new Tuple2<>(Long.parseLong(FilenameUtils.getBaseName(pair._1()).replace(SVConstants.FASTQ_OUT_PREFIX, "")), pair._2())) // essentially extract asm id from FASTQ file name
                                                                                .filter(pair -> interestingAsmIDs.contains(pair._1())) // filter out unused assemblies
                                                                                .mapValues(contents -> SVFastqUtils.convertToReads(SVFastqUtils.extractFASTQContents(contents))); // FASTQ-> GATKRead

        // TODO: add test for this complicated streaming
        // expand the pair rdd by expanding the assemblies to their associated junctions; also filters the reads
        final JavaPairRDD<SVJunction, List<GATKRead>> junctionReads = preFilteredReads.flatMapToPair( pair -> {
            final Long asmId = pair._1();
            final List<GATKRead> reads = pair._2();
            return asmid2SVJunctions.get(asmId).stream() // a list of junctions
                    .map( junc -> new Tuple2<>(junc, reads.stream().filter(read -> readSuitableForGenotypingJunction(junc, read)).collect(Collectors.toList()))) // each map to a pair(junction, filteredReads)
                    .collect(Collectors.toList()); // finally collect the pairs as a list
        });

        return junctionReads
                .groupByKey()
                .mapValues(it -> StreamSupport.stream(it.spliterator(), false).flatMap(List::stream).distinct().collect(Collectors.toList()) ); // collapse into a single list
    }

    @VisibleForTesting
    Map<Long, List<SVJunction>> getAssemblyToJunctionsMapping(final List<SVJunction> svJunctions) {
        // get mapping from assembly id to SVJunction's; see discussion on StackOverflow 38471056
        return svJunctions.stream()
                    .flatMap(svJunction -> svJunction.getAssemblyIDs().stream().map(asid -> new AbstractMap.SimpleEntry<>(asid, svJunction)))
                    .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(Map.Entry::getValue, Collectors.toList())));
    }

    @VisibleForTesting
    static List<GATKRead> mergeReadsByName(final List<GATKRead> fastqReads,
                                           final List<GATKRead> bamReads){
        // for reads appear in both, keep the one from bam
        final Set<String> common = Sets.intersection(fastqReads.stream().map(GATKRead::getName).collect(Collectors.toSet()),
                                                     bamReads.stream().map(GATKRead::getName).collect(Collectors.toSet()));
        bamReads.addAll( fastqReads.stream().filter(r -> !common.contains(r.getName())).collect(Collectors.toList()) );
        return bamReads;
    }

    @VisibleForTesting
    static List<VariantContext> collectAndSortGenotypedVC(final JavaPairRDD<SVJunction, GenotypeLikelihoods> genotypedJunctions,
                                                          final SAMSequenceDictionary referenceSequenceDictionary){
        return StreamSupport.stream(genotypedJunctions.map(pair -> pair._1().setPL(pair._2().getAsPLs())).collect().spliterator(), false)// set PL
                .map(SVJunction::getGenotypedVC) // add annotation and update VC
                .sorted((VariantContext v1, VariantContext v2) -> IntervalUtils.compareLocatables(v1, v2, referenceSequenceDictionary))
                .collect(Collectors.toList());
    }

    /**
     * TODO: should we emit reference confidence, we need more
     * Writes out new VCF FORMAT fields computed by this tool.
     */
    @VisibleForTesting
    static void output(final List<VariantContext> vcList,
                       final SAMSequenceDictionary referenceSequenceDictionary,
                       final GCSOptions options,
                       final String outFileName){

        try (final OutputStream outputStream = new BufferedOutputStream(BucketUtils.createFile(outFileName, options))){

            final VariantContextWriter vcfWriter = new VariantContextWriterBuilder()
                    .clearOptions()
                    .setOutputStream(outputStream)
                    .setReferenceDictionary(referenceSequenceDictionary)
                    .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .setOptions(Arrays.stream(new Options[]{}).collect(Collectors.toCollection(()->EnumSet.noneOf(Options.class))))
                    .setOption(Options.WRITE_FULL_FORMAT_FIELD)
                    .build();

            final VCFHeader header = new VCFHeader(GATKSVVCFHeaderLines.vcfHeaderLines.values().stream().collect(Collectors.toSet()), Collections.singletonList(testSampleName));
            header.setSequenceDictionary(referenceSequenceDictionary);
            header.addMetaDataLine( VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY) );

            header.addMetaDataLine( VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY) );
            header.addMetaDataLine( VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY));
            header.addMetaDataLine( VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_PL_KEY) );

            vcfWriter.writeHeader(header);
            vcList.forEach(vc -> logger.debug(String.format(".........%s:%d-%d\t%s.........", vc.getContig(), vc.getStart(), vc.getEnd(), vc.getGenotype(0).getLikelihoodsString())));
            vcList.forEach(vcfWriter::add);
            vcfWriter.close();
        } catch (final IOException ioex) {
            throw new GATKException("Could not create output file", ioex);
        }
    }

    // -----------------------------------------------------------------------------------------------
    // Utilities: debugging
    // -----------------------------------------------------------------------------------------------

    private void writeoutDebugReadLikelihoods(final JavaPairRDD<SVJunction, ReadLikelihoods<SVDummyAllele>> readLikelihoodsJavaPairRDD,
                                              final String debugStringSavingPath){
        readLikelihoodsJavaPairRDD.keys().map(j -> j.debugString).saveAsTextFile(debugStringSavingPath);
    }

}
