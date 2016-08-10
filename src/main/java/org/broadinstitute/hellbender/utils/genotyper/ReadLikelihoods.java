package org.broadinstitute.hellbender.utils.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.apache.commons.collections.ListUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Read-likelihoods container implementation based on integer indexed arrays.
 *
 * @param <A> the type of the allele the likelihood makes reference to.
 *
 * Note: this class uses FastUtil collections for speed.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadLikelihoods<A extends Allele> implements SampleList, AlleleList<A> {

    // -----------------------------------------------------------------------------------------------
    // Core structs
    // -----------------------------------------------------------------------------------------------

    /**
     * Indexed per sample, allele and finally read (within sample).
     * <p>
     *     valuesBySampleIndex[s][a][r] == lnLk(R_r | A_a) where R_r comes from Sample s.
     * </p>
     */
    private final double[][][] valuesBySampleIndex;

    /**
     * Contains information about the best allele for a read search result.
     */
    public final class BestAllele {
        public static final double INFORMATIVE_THRESHOLD = 0.2;

        /**
         * Null if there is no possible match (no allele?).
         */
        public final A allele;

        /**
         * The containing sample.
         */
        public final String sample;

        /**
         * The query read.
         */
        public final GATKRead read;

        /**
         * If {@link #allele} != null, then indicates the likelihood of the read.
         */
        public final double likelihood;

        /**
         * Confidence that the read actually was generated under that likelihood.
         * This is equal to the difference between this and the second best allele match.
         */
        public final double confidence;

        private BestAllele(final int sampleIndex, final int readIndex, final int bestAlleleIndex,
                           final double likelihood, final double secondBestLikelihood) {
            allele = bestAlleleIndex == -1 ? null : alleles.getAllele(bestAlleleIndex);
            sample = samples.getSample(sampleIndex);
            read = readsBySampleIndex[sampleIndex][readIndex];
            this.likelihood = likelihood;
            if(likelihood==Double.NEGATIVE_INFINITY){
                confidence = 0;
            }else{
                confidence = likelihood == secondBestLikelihood ? 0 : likelihood - secondBestLikelihood;
            }
        }

        public boolean isInformative() {
            return confidence > INFORMATIVE_THRESHOLD;
        }
    }

    /**
     * Implements a likelihood matrix per sample given its index.
     */
    private final class SampleMatrix implements LikelihoodMatrix<A> {

        private final int sampleIndex;

        private SampleMatrix(final int sampleIndex) {
            this.sampleIndex = sampleIndex;
        }

        @Override
        public List<GATKRead> reads() {
            return sampleReads(sampleIndex);
        }

        @Override
        public List<A> alleles() {
            return ReadLikelihoods.this.alleles();
        }

        @Override
        public void set(final int alleleIndex, final int readIndex, final double value) {
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            Utils.validIndex(readIndex, valuesBySampleIndex[sampleIndex][alleleIndex].length);
            valuesBySampleIndex[sampleIndex][alleleIndex][readIndex] = value;
        }

        @Override
        public double get(final int alleleIndex, final int readIndex) {
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            Utils.validIndex(readIndex, valuesBySampleIndex[sampleIndex][alleleIndex].length);
            return valuesBySampleIndex[sampleIndex][alleleIndex][readIndex];
        }

        @Override
        public int indexOfAllele(final A allele) {
            Utils.nonNull(allele);
            return ReadLikelihoods.this.indexOfAllele(allele);
        }

        @Override
        public int indexOfRead(final GATKRead read) {
            Utils.nonNull(read);
            return ReadLikelihoods.this.readIndex(sampleIndex, read);
        }

        @Override
        public int numberOfAlleles() {
            return alleles.numberOfAlleles();
        }

        @Override
        public int numberOfReads() {
            return readsBySampleIndex[sampleIndex].length;
        }

        @Override
        public A getAllele(final int alleleIndex) {
            return ReadLikelihoods.this.getAllele(alleleIndex);
        }

        @Override
        public GATKRead getRead(final int readIndex) {
            final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
            Utils.validIndex(readIndex, sampleReads.length);
            return sampleReads[readIndex];
        }

        @Override
        public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
            Utils.nonNull(dest);
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            System.arraycopy(valuesBySampleIndex[sampleIndex][alleleIndex], 0, dest, offset, numberOfReads());
        }
    }

    /**
     * Sample matrices lazily initialized (the elements not the array) by invoking {@link #sampleMatrix(int)}.
     */
    private final LikelihoodMatrix<A>[] sampleMatrices;

    // -----------------------------------------------------------------------------------------------
    // Trivial structs
    // -----------------------------------------------------------------------------------------------

    /**
     * Reads by sample index. Each sub array contains reference to the reads of the ith sample.
     */
    private final GATKRead[][] readsBySampleIndex;

    /**
     * Sample list
     */
    private final SampleList samples;

    /**
     * Allele list
     */
    private AlleleList<A> alleles;

    /**
     * Cached allele list.
     */
    private List<A> alleleList;

    /**
     * Cached sample list.
     */
    private List<String> sampleList;

    /**
     * Maps from each read to its index within the sample.
     *
     * <p>In order to save CPU time the indices contained in this array (not the array itself) is
     * lazily initialized by invoking {@link #readIndicesBySampleIndex(int)}.</p>
     */
    private final Object2IntMap<GATKRead>[] readIndexBySampleIndex;

    /**
     * Index indicating that the reference allele is missing.
     */
    private static final int MISSING_REF = -1;

    /**
     * Index of the reference allele if any, otherwise {@link #MISSING_REF}.
     */
    private int referenceAlleleIndex = MISSING_REF;

    /**
     * Caches the read-list per sample list returned by {@link #sampleReads(int)}
     */
    private final List<GATKRead>[] readListBySampleIndex;

    // -----------------------------------------------------------------------------------------------
    // CORE (balanced modifiers)
    // -----------------------------------------------------------------------------------------------

    /**
     * Adjusts likelihoods so that for each read, the best allele likelihood is 0 and caps the minimum likelihood
     * of any allele for each read based on the maximum alternative allele likelihood.
     *
     * @param bestToZero                        set the best likelihood to 0, others will be subtracted the same amount.
     * @param maximumLikelihoodDifferenceCap    maximum difference between the best alternative allele likelihood
     *                                           and any other likelihood.
     *
     * @throws IllegalArgumentException if {@code maximumDifferenceWithBestAlternative} is not 0 or less.
     */
    public void normalizeLikelihoods(final boolean bestToZero,
                                     final double maximumLikelihoodDifferenceCap) {
        Utils.validateArg(maximumLikelihoodDifferenceCap < 0.0 && !Double.isNaN(maximumLikelihoodDifferenceCap),
                "the minimum reference likelihood fall must be negative");

        if (maximumLikelihoodDifferenceCap == Double.NEGATIVE_INFINITY && !bestToZero) {
            return;
        }

        final int alleleCount = alleles.numberOfAlleles();
        if (alleleCount == 0){ // trivial case there is no alleles.
            return;
        } else if (alleleCount == 1 && !bestToZero) {
            return;
        }

        for (int s = 0; s < valuesBySampleIndex.length; s++) {
            final double[][] sampleValues = valuesBySampleIndex[s];
            final int readCount = readsBySampleIndex[s].length;
            for (int r = 0; r < readCount; r++) {
                normalizeLikelihoodsPerRead(bestToZero, maximumLikelihoodDifferenceCap, sampleValues, s, r);
            }
        }
    }

    /**
     * Does the normalizeLikelihoods job for each read.
     * @param bestToZero                            set the best likelihood to 0, others will be subtracted the same amount.
     * @param maximumBestAltLikelihoodDifference    maximum difference between the best alternative allele likelihood and any other likelihood.
     * @param sampleValues                          loglikelihood matrix of a sample (indexed by {@code sampleIndex}
     * @param sampleIndex                           index of sample
     * @param readIndex                             index of read into the particular sample
     */
    private void normalizeLikelihoodsPerRead(final boolean bestToZero,
                                             final double maximumBestAltLikelihoodDifference,
                                             final double[][] sampleValues,
                                             final int sampleIndex,
                                             final int readIndex) {

        // first find the best alt allele
        final BestAllele bestAlternativeAllele = searchBestAllele(sampleIndex,readIndex,false);

        // left tail
        final double worstLikelihoodCap = bestAlternativeAllele.likelihood + maximumBestAltLikelihoodDifference;

        // find reference allele likelihood
        final double referenceLikelihood = referenceAlleleIndex == MISSING_REF ? Double.NEGATIVE_INFINITY :
                                                                                 sampleValues[referenceAlleleIndex][readIndex];

        // get the max out of ref allele likelihood and best alt allele likelihood
        final double bestAbsoluteLikelihood = Math.max(bestAlternativeAllele.likelihood, referenceLikelihood);

        final int alleleCount = alleles.numberOfAlleles();
        if (bestToZero) {
            if (bestAbsoluteLikelihood == Double.NEGATIVE_INFINITY) { // if all alleles are bad, shift all likelihood to zero
                for (int a = 0; a < alleleCount; a++) {
                    sampleValues[a][readIndex] = 0;
                }
            } else if (worstLikelihoodCap != Double.NEGATIVE_INFINITY) {
                for (int a = 0; a < alleleCount; a++) { // first cap all allele's likelihood from below (worstLikelihoodcap) then shift down by best likelihood -> max would become 0
                    sampleValues[a][readIndex] = (sampleValues[a][readIndex] < worstLikelihoodCap ? worstLikelihoodCap : sampleValues[a][readIndex]) - bestAbsoluteLikelihood;
                }
            } else { // shift down all likelihoods when ref is the best and all alt allele are bad (best alt allele has -\infty likelihood)
                for (int a = 0; a < alleleCount; a++) {
                    sampleValues[a][readIndex] -= bestAbsoluteLikelihood;
                }
            }
        } else {
            // Guarantee to be the case by enclosing code.
            for (int a = 0; a < alleleCount; a++) {
                if (sampleValues[a][readIndex] < worstLikelihoodCap) { // caps all alleles' (including ref if it's not the best) likelihood from below
                    sampleValues[a][readIndex] = worstLikelihoodCap;
                }
            }
        }
    }

    /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each read in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @param overlap           if not {@code null}, only reads that overlap the location (with unclipping) will be present
     *                          in the output read-collection.
     *
     * @return                  never {@code null}. The result will have the requested set of new alleles
     *                          (keys in {@code newToOldAlleleMap}, and the same set of samples and reads as the original.
     *
     * @throws IllegalArgumentException if
     *                                  1) {@code newToOldAlleleMap} is {@code null}, or
     *                                  2) contains {@code null} values, or
     *                                  3) its values contain reference to non-existing alleles in this read-likelihood collection.
     *                                  4) Also no new allele can have zero old alleles mapping
     *                                     nor two new alleles can make reference to the same old allele.
     */
    public <B extends Allele> ReadLikelihoods<B> marginalize(final Map<B, List<A>> newToOldAlleleMap,
                                                             final Locatable overlap) {
        Utils.nonNull(newToOldAlleleMap, "the input allele mapping cannot be null");
        if (overlap == null) {
            return marginalize(newToOldAlleleMap);
        }

        @SuppressWarnings("unchecked")
        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final int[] oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        final int[][] readsToKeep = overlappingReadIndicesBySampleIndex(overlap);

        // We calculate the marginal likelihoods.
        final double[][][] newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, oldToNewAlleleIndexMap, readsToKeep);

        final int sampleCount = samples.numberOfSamples();

        @SuppressWarnings({"rawtypes","unchecked"})
        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = (Object2IntMap<GATKRead>[])new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];

        for (int s = 0; s < sampleCount; s++) {
            final int[] sampleReadsToKeep = readsToKeep[s];
            final GATKRead[] oldSampleReads = readsBySampleIndex[s];
            final int oldSampleReadCount = oldSampleReads.length;
            final int newSampleReadCount = sampleReadsToKeep.length;
            if (newSampleReadCount == oldSampleReadCount) {
                newReadsBySampleIndex[s] = oldSampleReads.clone();
            } else {
                newReadsBySampleIndex[s] = new GATKRead[newSampleReadCount];
                for (int i = 0; i < newSampleReadCount; i++) {
                    newReadsBySampleIndex[s][i] = oldSampleReads[sampleReadsToKeep[i]];
                }
            }
        }

        // Finally we create the new read-likelihood
        return new ReadLikelihoods<>(new IndexedAlleleList<>(newAlleles), samples, newReadsBySampleIndex, newReadIndexBySampleIndex, newLikelihoodValues);
    }

    /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each read in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @return never {@code null}. The result will have the requested set of new alleles (keys in {@code newToOldAlleleMap}, and
     * the same set of samples and reads as the original.
     *
     * @throws IllegalArgumentException is {@code newToOldAlleleMap} is {@code null} or contains {@code null} values,
     *  or its values contain reference to non-existing alleles in this read-likelihood collection. Also no new allele
     *  can have zero old alleles mapping nor two new alleles can make reference to the same old allele.
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <B extends Allele> ReadLikelihoods<B> marginalize(final Map<B, List<A>> newToOldAlleleMap) {
        Utils.nonNull(newToOldAlleleMap, "the input allele mapping cannot be null");

        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final int[] oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        // We calculate the marginal likelihoods.
        final double[][][] newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, oldToNewAlleleIndexMap, null);

        final int sampleCount = samples.numberOfSamples();

        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];

        for (int s = 0; s < sampleCount; s++) {
            newReadsBySampleIndex[s] = readsBySampleIndex[s].clone();
        }

        // Finally we create the new read-likelihood
        return new ReadLikelihoods<>(new IndexedAlleleList<>(newAlleles), samples, newReadsBySampleIndex, newReadIndexBySampleIndex, newLikelihoodValues);
    }

    /**
     * Calculate the marginal likelihoods considering the old -> new allele index mapping.
     * @param oldAlleleCount
     * @param newAlleleCount
     * @param oldToNewAlleleIndexMap
     * @param readsToKeep
     * @return
     */
    private double[][][] marginalLikelihoods(final int oldAlleleCount,
                                             final int newAlleleCount,
                                             final int[] oldToNewAlleleIndexMap,
                                             final int[][] readsToKeep) {

        final int sampleCount = samples.numberOfSamples();
        final double[][][] result = new double[sampleCount][][];

        for (int s = 0; s < sampleCount; s++) {

            final int oldSampleReadCount = readsBySampleIndex[s].length;
            final int[] sampleReadToKeep = readsToKeep == null || readsToKeep[s].length == oldSampleReadCount ? null : readsToKeep[s];
            final int newSampleReadCount = sampleReadToKeep == null ? oldSampleReadCount : sampleReadToKeep.length;

            // old result and new target initialization
            final double[][] oldSampleValues = valuesBySampleIndex[s];
            final double[][] newSampleValues = result[s] = new double[newAlleleCount][newSampleReadCount];
            // We initiate all likelihoods to -Inf.
            for (int a = 0; a < newAlleleCount; a++) {
                Arrays.fill(newSampleValues[a], Double.NEGATIVE_INFINITY);
            }

            // For each old allele and read we update the new table keeping the maximum likelihood.
            for (int r = 0; r < newSampleReadCount; r++) {
                for (int a = 0; a < oldAlleleCount; a++) {
                    final int oldReadIndex = newSampleReadCount == oldSampleReadCount ? r : sampleReadToKeep[r];
                    final int newAlleleIndex = oldToNewAlleleIndexMap[a];
                    if (newAlleleIndex == -1) {
                        continue;
                    }
                    final double likelihood = oldSampleValues[a][oldReadIndex];
                    if (likelihood > newSampleValues[newAlleleIndex][r]) {
                        newSampleValues[newAlleleIndex][r] = likelihood;
                    }
                }
            }
        }

        return result;
    }

    /**
     * Change reads based on input mapping from old read to new read.
     * @param readRealignments
     */
    public void changeReads(final Map<GATKRead, GATKRead> readRealignments) {
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final Object2IntMap<GATKRead> readIndex = readIndexBySampleIndex[s];
            final int sampleReadCount = sampleReads.length;
            for (int r = 0; r < sampleReadCount; r++) {
                final GATKRead read = sampleReads[r];
                final GATKRead replacement = readRealignments.get(read);
                if (replacement == null) {
                    continue;
                }
                sampleReads[r] = replacement;
                if (readIndex != null) {
                    readIndex.remove(read);
                    readIndex.put(replacement, r);
                }
            }
        }
    }

    // -----------------------------------------------------------------------------------------------
    // Serving core
    // -----------------------------------------------------------------------------------------------

    /**
     * Returns a 2-D array with row-dimension representing samples, i.e. row count will be the same as sample count.
     * The sub arrays will be a sample's reads' index that overlaps (with unclipping) the provided interval.
     *
     * @return {@code null} if the provided {@code overlap} is null
     */
    private int[][] overlappingReadIndicesBySampleIndex(final Locatable overlap) {
        if (overlap == null) {
            return null;
        }
        final int sampleCount = samples.numberOfSamples();
        final int[][] result = new int[sampleCount][];
        final IntArrayList buffer = new IntArrayList(200);
        for (int s = 0; s < sampleCount; s++) {
            buffer.clear();
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final int sampleReadCount = sampleReads.length;
            buffer.ensureCapacity(sampleReadCount);
            for (int r = 0; r < sampleReadCount; r++) {
                if (unclippedReadOverlapsRegion(sampleReads[r], overlap)) {
                    buffer.add(r);
                }
            }
            result[s] = buffer.toIntArray();
        }
        return result;
    }

    /**
     * calculates an old to new allele index map array.
     * -1 indicates the old allele doesn't map to any new allele.
     * @param newToOldAlleleMap
     * @param oldAlleleCount
     * @param newAlleles
     * @param <B>
     * @return
     */
    private <B extends Allele> int[] oldToNewAlleleIndexMap(final Map<B, List<A>> newToOldAlleleMap,
                                                            final int oldAlleleCount,
                                                            final B[] newAlleles) {
        Arrays.stream(newAlleles).forEach(Utils::nonNull);
        Utils.containsNoNull(newToOldAlleleMap.values(), "no new allele list can be null");
        newToOldAlleleMap.values().stream().forEach(oldList -> Utils.containsNoNull(oldList,"old alleles cannot be null"));

        final int[] oldToNewAlleleIndexMap = new int[oldAlleleCount];
        Arrays.fill(oldToNewAlleleIndexMap, -1);  // -1 indicate that there is no new allele that make reference to that old one.

        for (int newIndex = 0; newIndex < newAlleles.length; newIndex++) {
            final B newAllele = newAlleles[newIndex];
            for (final A oldAllele : newToOldAlleleMap.get(newAllele)) {
                final int oldAlleleIndex = indexOfAllele(oldAllele);
                if (oldAlleleIndex == -1) {
                    throw new IllegalArgumentException("missing old allele " + oldAllele + " in likelihood collection ");
                }
                if (oldToNewAlleleIndexMap[oldAlleleIndex] != -1) {
                    throw new IllegalArgumentException("collision: two new alleles make reference to the same old allele");
                }
                oldToNewAlleleIndexMap[oldAlleleIndex] = newIndex;
            }
        }
        return oldToNewAlleleIndexMap;
    }

    // -----------------------------------------------------------------------------------------------
    // Initializers
    // -----------------------------------------------------------------------------------------------

    /**
     * Constructs a new read-likelihood collection.
     *
     * <p>
     *     The initial likelihoods for all allele-read combinations are
     *     0.
     * </p>
     *
     * @param samples all supported samples in the collection.
     * @param alleles all supported alleles in the collection.
     * @param reads reads stratified per sample.
     *
     * @throws IllegalArgumentException if any of {@code allele}, {@code samples}
     * or {@code reads} is {@code null},
     *  or if they contain null values.
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    public ReadLikelihoods(final SampleList samples,
                           final AlleleList<A> alleles,
                           final Map<String, List<GATKRead>> reads) {
        Utils.nonNull(alleles, "allele list cannot be null");
        Utils.nonNull(samples, "sample list cannot be null");
        Utils.nonNull(reads, "read map cannot be null");

        this.samples = samples;
        this.alleles = alleles;

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        readsBySampleIndex = new GATKRead[sampleCount][];
        readListBySampleIndex = (List<GATKRead>[])new List[sampleCount];
        valuesBySampleIndex = new double[sampleCount][][];
        referenceAlleleIndex = findReferenceAllele(alleles);

        readIndexBySampleIndex = new Object2IntMap[sampleCount];

        setupIndexes(reads, sampleCount, alleleCount);

        sampleMatrices = (LikelihoodMatrix<A>[]) new LikelihoodMatrix[sampleCount];
    }

    // Internally used constructor.
    @SuppressWarnings({"unchecked", "rawtypes"})
    private ReadLikelihoods(final AlleleList alleles,
                            final SampleList samples,
                            final GATKRead[][] readsBySampleIndex,
                            final Object2IntMap<GATKRead>[] readIndex,
                            final double[][][] values) {
        this.samples = samples;
        this.alleles = alleles;
        this.readsBySampleIndex = readsBySampleIndex;
        this.valuesBySampleIndex = values;
        this.readIndexBySampleIndex = readIndex;
        final int sampleCount = samples.numberOfSamples();
        this.readListBySampleIndex = (List<GATKRead>[])new List[sampleCount];

        referenceAlleleIndex = findReferenceAllele(alleles);
        sampleMatrices = (LikelihoodMatrix<A>[]) new LikelihoodMatrix[sampleCount];
    }

    /**
     *  Add all the indices to alleles, sample and reads in the look-up maps.
     * @param reads
     * @param sampleCount
     * @param alleleCount
     */
    void setupIndexes(final Map<String, List<GATKRead>> reads,
                      final int sampleCount,
                      final int alleleCount) {
        for (int i = 0; i < sampleCount; i++) {
            setupSampleData(i, reads, alleleCount);
        }
    }

    /**
     * Assumes that {@link #samples} has been initialized with the sample names.
     * @param sampleIndex
     * @param readsBySample
     * @param alleleCount
     */
    private void setupSampleData(final int sampleIndex,
                                 final Map<String, List<GATKRead>> readsBySample,
                                 final int alleleCount) {
        final String sample = samples.getSample(sampleIndex);

        final List<GATKRead> reads = readsBySample.get(sample);
        readsBySampleIndex[sampleIndex] = reads == null
                ? new GATKRead[0]
                : reads.toArray(new GATKRead[reads.size()]);
        final int sampleReadCount = readsBySampleIndex[sampleIndex].length;

        final double[][] sampleValues = new double[alleleCount][sampleReadCount];
        valuesBySampleIndex[sampleIndex] = sampleValues;
    }

    /**
     * Search for the reference allele in the provided list of alleles.
     * @param alleles   a list of alleles in which the search is to be performed
     * @return          index of the reference allele in the input;
     *                  {@link #MISSING_REF} if not found
     */
    private static int findReferenceAllele(final AlleleList<?> alleles) {
        return IntStream.range(0, alleles.numberOfAlleles()).filter(i -> alleles.getAllele(i).isReference()).findAny().orElse(MISSING_REF);
    }

    // -----------------------------------------------------------------------------------------------
    // Modifiers (modifies fields) : add
    // -----------------------------------------------------------------------------------------------

    /**
     * Adds the non-reference allele to the read-likelihood collection setting each read likelihood to the second
     * best found (or best one if only one allele has likelihood).
     *
     * <p>Nothing will happen if the read-likelihoods collection already includes the non-ref allele</p>
     *
     * <p>
     *     <i>Implementation note: even when strictly speaking we do not need to demand the calling code to pass
     *     the reference the non-ref allele, we still demand it in order to lead the
     *     the calling code to use the right generic type for this likelihoods
     *     collection {@link Allele}.</i>
     * </p>
     *
     * @param nonRefAllele the non-ref allele.
     *
     * @throws IllegalArgumentException if {@code nonRefAllele} is anything but the designated &lt;NON_REF&gt;
     * symbolic allele {@link org.broadinstitute.hellbender.utils.variant.GATKVCFConstants#NON_REF_SYMBOLIC_ALLELE}.
     */
    public void addNonReferenceAllele(final A nonRefAllele) {
        Utils.nonNull(nonRefAllele, "non-ref allele cannot be null");
        if (!nonRefAllele.equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)) {
            throw new IllegalArgumentException("the non-ref allele is not valid");
        }
        if (alleles.containsAllele(nonRefAllele)) {
            return;
        }

        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = oldAlleleCount + 1;
        @SuppressWarnings("unchecked")
        final A[] newAlleles = (A[]) new Allele[newAlleleCount];
        for (int a = 0; a < oldAlleleCount; a++) {
            newAlleles[a] = alleles.getAllele(a);
        }
        newAlleles[oldAlleleCount] = nonRefAllele;
        alleles = new IndexedAlleleList<>(newAlleles);
        alleleList = null; // remove the cached alleleList.

        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            addNonReferenceAlleleLikelihoodsPerSample(oldAlleleCount, newAlleleCount, s);
        }
    }

    /**
     * Add alleles that are missing in the read-likelihoods collection giving all reads a default
     * likelihood value.
     * @param candidateAlleles the potentially missing alleles.
     * @param defaultLikelihood the default read likelihood value for that allele.
     *
     * @throws IllegalArgumentException if {@code candidateAlleles} is {@code null} or there is more than
     * one missing allele that is a reference or there is one but the collection already has
     * a reference allele.
     */
    public void addMissingAlleles(final Collection<A> candidateAlleles,
                                  final double defaultLikelihood) {
        Utils.nonNull(candidateAlleles, "the candidateAlleles list cannot be null");
        if (candidateAlleles.isEmpty()) {
            return;
        }
        final List<A> allelesToAdd = candidateAlleles.stream().filter(allele -> !alleles.containsAllele(allele)).collect(Collectors.toList());

        if (allelesToAdd.isEmpty()) {
            return;
        }

        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = alleles.numberOfAlleles() + allelesToAdd.size();

        alleleList = null;
        int referenceIndex = this.referenceAlleleIndex;

        @SuppressWarnings("unchecked")
        final List<A> newAlleles = ListUtils.union(alleles.asListOfAlleles(), allelesToAdd);
        alleles = new IndexedAlleleList<>(newAlleles);

        // if we previously had no reference allele, update the reference index if a reference allele is added
        // if we previously had a reference and try to add another, throw an exception
        final OptionalInt indexOfReferenceInAllelesToAdd = IntStream.range(0, allelesToAdd.size())
                .filter(n -> allelesToAdd.get(n).isReference()).findFirst();
        if (referenceIndex != MISSING_REF) {
            Utils.validateArg(!indexOfReferenceInAllelesToAdd.isPresent(), "there can only be one reference allele");
        } else if (indexOfReferenceInAllelesToAdd.isPresent()){
            referenceAlleleIndex = oldAlleleCount + indexOfReferenceInAllelesToAdd.getAsInt();
        }

        //copy old allele likelihoods and set new allele likelihoods to the default value
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final int sampleReadCount = readsBySampleIndex[s].length;
            final double[][] newValuesBySampleIndex = Arrays.copyOf(valuesBySampleIndex[s], newAlleleCount);
            for (int a = oldAlleleCount; a < newAlleleCount; a++) {
                newValuesBySampleIndex[a] = new double[sampleReadCount];
                if (defaultLikelihood != 0.0) {
                    Arrays.fill(newValuesBySampleIndex[a], defaultLikelihood);
                }
            }
            valuesBySampleIndex[s] = newValuesBySampleIndex;
        }
    }

    /**
     * Add more reads to the collection.
     *
     * @param readsBySample reads to add.
     * @param initialLikelihood the likelihood for the new entries.
     *
     * @throws IllegalArgumentException if {@code readsBySample} is {@code null} or {@code readsBySample} contains
     *  {@code null} reads, or {@code readsBySample} contains read that are already present in the read-likelihood
     *  collection.
     */
    public void addReads(final Map<String,List<GATKRead>> readsBySample,
                         final double initialLikelihood) {
        for (final Map.Entry<String,List<GATKRead>> entry : readsBySample.entrySet()) {
            final String sample = entry.getKey();
            final List<GATKRead> newSampleReads = entry.getValue();
            final int sampleIndex = samples.indexOfSample(sample);

            if (sampleIndex == -1) {
                throw new IllegalArgumentException("input sample " + sample +
                        " is not part of the read-likelihoods collection");
            }

            if (newSampleReads == null || newSampleReads.isEmpty()) {
                continue;
            }

            final int sampleReadCount = readsBySampleIndex[sampleIndex].length;
            final int newSampleReadCount = sampleReadCount + newSampleReads.size();

            appendReads(newSampleReads, sampleIndex, sampleReadCount, newSampleReadCount);
            extendsLikelihoodArrays(initialLikelihood, sampleIndex, sampleReadCount, newSampleReadCount);
        }
    }


    /**
     * Append the new read reference into the structure per-sample.
     * @param newSampleReads
     * @param sampleIndex
     * @param sampleReadCount
     * @param newSampleReadCount
     */
    private void appendReads(final List<GATKRead> newSampleReads,
                             final int sampleIndex,
                             final int sampleReadCount,
                             final int newSampleReadCount) {
        final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex] =
                Arrays.copyOf(readsBySampleIndex[sampleIndex], newSampleReadCount);

        int nextReadIndex = sampleReadCount;
        final Object2IntMap<GATKRead> sampleReadIndex = readIndexBySampleIndex[sampleIndex];
        for (final GATKRead newRead : newSampleReads) {
            //    if (sampleReadIndex.containsKey(newRead)) // might be worth handle this without exception (ignore the read?) but in practice should never be the case.
            //        throw new IllegalArgumentException("you cannot add reads that are already in read-likelihood collection");
            if (sampleReadIndex != null ) {
                sampleReadIndex.put(newRead, nextReadIndex);
            }
            sampleReads[nextReadIndex++] = newRead;
        }
    }

    /**
     * Extends the likelihood arrays-matrices.
     * @param initialLikelihood
     * @param sampleIndex
     * @param sampleReadCount
     * @param newSampleReadCount
     */
    private void extendsLikelihoodArrays(final double initialLikelihood,
                                         final int sampleIndex,
                                         final int sampleReadCount,
                                         final int newSampleReadCount) {
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        final int alleleCount = alleles.numberOfAlleles();
        for (int a = 0; a < alleleCount; a++) {
            sampleValues[a] = Arrays.copyOf(sampleValues[a], newSampleReadCount);
        }
        if (initialLikelihood != 0.0) // the default array new value.
        {
            for (int a = 0; a < alleleCount; a++) {
                Arrays.fill(sampleValues[a], sampleReadCount, newSampleReadCount, initialLikelihood);
            }
        }
    }

    /**
     * Updates per-sample structures according to the addition of the NON_REF allele.
     * @param alleleCount
     * @param newAlleleCount
     * @param sampleIndex
     */
    private void addNonReferenceAlleleLikelihoodsPerSample(final int alleleCount,
                                                           final int newAlleleCount,
                                                           final int sampleIndex) {
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex] = Arrays.copyOf(valuesBySampleIndex[sampleIndex], newAlleleCount);
        final int sampleReadCount = readsBySampleIndex[sampleIndex].length;

        final double[] nonRefAlleleLikelihoods = sampleValues[alleleCount] = new double [sampleReadCount];
        Arrays.fill(nonRefAlleleLikelihoods, Double.NEGATIVE_INFINITY);
        for (int r = 0; r < sampleReadCount; r++) {
            final BestAllele bestAllele = searchBestAllele(sampleIndex,r,true);
            final double secondBestLikelihood = Double.isInfinite(bestAllele.confidence) ? bestAllele.likelihood
                    : bestAllele.likelihood - bestAllele.confidence;
            nonRefAlleleLikelihoods[r] = secondBestLikelihood;
        }
    }

    // -----------------------------------------------------------------------------------------------
    // Modifiers (modifies fields) : remove
    // -----------------------------------------------------------------------------------------------

    /**
     * Removes those read that the best possible likelihood given any allele is just too low.
     *
     * <p>
     *     This is determined by a maximum error per read-base against the best likelihood possible.
     * </p>
     *
     * @param maximumErrorPerBase the minimum acceptable error rate per read base, must be
     *                            a positive number.
     *
     * @throws IllegalStateException is not supported for read-likelihood that do not contain alleles.
     *
     * @throws IllegalArgumentException if {@code maximumErrorPerBase} is negative.
     */
    public void filterPoorlyModeledReads(final double maximumErrorPerBase) {
        Utils.validateArg(alleles.numberOfAlleles() > 0, "unsupported for read-likelihood collections with no alleles");
        Utils.validateArg(!Double.isNaN(maximumErrorPerBase) && maximumErrorPerBase > 0.0, "the maximum error per base must be a positive number");

        new IndexRange(0, samples.numberOfSamples()).forEach(s -> {
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final List<Integer> removeIndices = new IndexRange(0, sampleReads.length).filter(r -> readIsPoorlyModelled(s, r, sampleReads[r], maximumErrorPerBase));
            removeSampleReads(s, removeIndices, alleles.numberOfAlleles());
        });

    }

    /**
     * Model is the same as:
     * @see {@link PerReadAlleleLikelihoodMap#readIsPoorlyModelled(GATKRead, Collection, double)}.
     */
    private boolean readIsPoorlyModelled(final int sampleIndex,
                                         final int readIndex,
                                         final GATKRead read,
                                         final double maxErrorRatePerBase) {
        final double maxErrorsForRead = Math.min(2.0, Math.ceil(read.getLength() * maxErrorRatePerBase));
        final double log10QualPerBase = -4.0;
        final double log10MaxLikelihoodForTrueAllele = maxErrorsForRead * log10QualPerBase;

        final int alleleCount = alleles.numberOfAlleles();
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][readIndex] >= log10MaxLikelihoodForTrueAllele) {
                return false;
            }
        }
        return true;
    }

    /**
     * Remove those reads that do not overlap certain genomic location.
     *
     * <p>
     *     This method modifies the current read-likelihoods collection.
     * </p>
     *
     * @param location the target location.
     *
     * @throws IllegalArgumentException the location cannot be {@code null} nor unmapped.
     */
    public void filterToOnlyOverlappingUnclippedReads(final SimpleInterval location) {
        Utils.nonNull(location, "the location cannot be null");

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();
        for (int s = 0; s < sampleCount; s++) {
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final List<Integer> removeIndices = new IndexRange(0, sampleReads.length)
                    .filter(r -> !unclippedReadOverlapsRegion(sampleReads[r], location));
            removeSampleReads(s, removeIndices, alleleCount);
        }
    }

    /**
     * Downsamples reads based on contamination fractions making sure that all alleles are affected proportionally.
     *
     * @param perSampleDownsamplingFraction contamination sample map where the sample name are the keys and the
     *                                       fractions are the values.
     *
     * @throws IllegalArgumentException if {@code perSampleDownsamplingFraction} is {@code null}.
     */
    public void contaminationDownsampling(final Map<String, Double> perSampleDownsamplingFraction) {
        Utils.nonNull(perSampleDownsamplingFraction);

        final int alleleCount = alleles.numberOfAlleles();
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final String sample = samples.getSample(s);
            final Double fractionDouble = perSampleDownsamplingFraction.get(sample);
            if (fractionDouble == null) {
                continue;
            }
            final double fraction = fractionDouble;
            if (Double.isNaN(fraction) || fraction <= 0.0) {
                continue;
            }
            if (fraction >= 1.0) {
                final List<Integer> removeIndices = IntStream.range(0, readsBySampleIndex[s].length).boxed().collect(Collectors.toList());
                removeSampleReads(s, removeIndices, alleleCount);
            } else {
                final Map<A,List<GATKRead>> readsByBestAllelesMap = readsByBestAlleleMap(s);
                removeSampleReads(s, AlleleBiasedDownsamplingUtils.selectAlleleBiasedReads(readsByBestAllelesMap, fraction),alleleCount);
            }
        }
    }

    private void removeSampleReads(final int sampleIndex,
                                   final List<Integer> removeIndices,
                                   final int alleleCount) {
        if (removeIndices.isEmpty()) {
            return;
        }

        final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
        final int sampleReadCount = sampleReads.length;

        final Object2IntMap<GATKRead> indexByRead = readIndexBySampleIndex[sampleIndex];
        if (indexByRead != null) {
            removeIndices.stream().forEach(n -> indexByRead.remove(sampleReads[n]));
        }
        final boolean[] removeIndex = new boolean[sampleReadCount];
        final int firstDeleted = removeIndices.get(0);
        removeIndices.stream().forEach(n -> removeIndex[n] = true);

        final int newSampleReadCount = sampleReadCount - removeIndices.size();

        // Now we skim out the removed reads from the read array.
        final GATKRead[] oldSampleReads = readsBySampleIndex[sampleIndex];
        final GATKRead[] newSampleReads = new GATKRead[newSampleReadCount];

        System.arraycopy(oldSampleReads, 0, newSampleReads, 0, firstDeleted);
        Utils.skimArray(oldSampleReads,firstDeleted, newSampleReads, firstDeleted, removeIndex, firstDeleted);

        // Then we skim out the likelihoods of the removed reads.
        final double[][] oldSampleValues = valuesBySampleIndex[sampleIndex];
        final double[][] newSampleValues = new double[alleleCount][newSampleReadCount];
        for (int a = 0; a < alleleCount; a++) {
            System.arraycopy(oldSampleValues[a], 0, newSampleValues[a], 0, firstDeleted);
            Utils.skimArray(oldSampleValues[a], firstDeleted, newSampleValues[a], firstDeleted, removeIndex, firstDeleted);
        }
        valuesBySampleIndex[sampleIndex] = newSampleValues;
        readsBySampleIndex[sampleIndex] = newSampleReads;
        readListBySampleIndex[sampleIndex] = null; // reset the unmodifiable list.
    }

    /**
     * Requires that the collection passed iterator can remove elements, and it can be modified.
     * @param sampleIndex
     * @param readsToRemove
     * @param alleleCount
     */
    private void removeSampleReads(final int sampleIndex,
                                   final Collection<GATKRead> readsToRemove,
                                   final int alleleCount) {
        final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
        final int sampleReadCount = sampleReads.length;

        final Object2IntMap<GATKRead> indexByRead = readIndicesBySampleIndex(sampleIndex);
        // Count how many we are going to remove, which ones (indexes) and remove entry from the read-index map.
        final boolean[] removeIndex = new boolean[sampleReadCount];
        int removeCount = 0; // captures the number of deletions.
        int firstDeleted = sampleReadCount;    // captures the first position that was deleted.

        final Iterator<GATKRead> readsToRemoveIterator = readsToRemove.iterator();
        while (readsToRemoveIterator.hasNext()) {
            final GATKRead read = readsToRemoveIterator.next();
            if (indexByRead.containsKey(read)) {
                final int index = indexByRead.getInt(read);
                if (firstDeleted > index) {
                    firstDeleted = index;
                }
                removeCount++;
                removeIndex[index] = true;
                readsToRemoveIterator.remove();
                indexByRead.remove(read);
            }
        }

        // Nothing to remove we just finish here.
        if (removeCount == 0) {
            return;
        }

        final int newSampleReadCount = sampleReadCount - removeCount;

        // Now we skim out the removed reads from the read array.
        final GATKRead[] oldSampleReads = readsBySampleIndex[sampleIndex];
        final GATKRead[] newSampleReads = new GATKRead[newSampleReadCount];

        System.arraycopy(oldSampleReads,0,newSampleReads,0,firstDeleted);
        Utils.skimArray(oldSampleReads,firstDeleted, newSampleReads, firstDeleted, removeIndex, firstDeleted);

        // Update the indices for the extant reads from the first deletion onwards.
        for (int r = firstDeleted; r < newSampleReadCount; r++) {
            indexByRead.put(newSampleReads[r], r);
        }

        // Then we skim out the likelihoods of the removed reads.
        final double[][] oldSampleValues = valuesBySampleIndex[sampleIndex];
        final double[][] newSampleValues = new double[alleleCount][newSampleReadCount];
        for (int a = 0; a < alleleCount; a++) {
            System.arraycopy(oldSampleValues[a],0,newSampleValues[a],0,firstDeleted);
            Utils.skimArray(oldSampleValues[a], firstDeleted, newSampleValues[a], firstDeleted, removeIndex, firstDeleted);
        }
        valuesBySampleIndex[sampleIndex] = newSampleValues;
        readsBySampleIndex[sampleIndex] = newSampleReads;
        readListBySampleIndex[sampleIndex] = null; // reset the unmodifiable list.
    }

    // -----------------------------------------------------------------------------------------------
    // Utilities (accessors for other methods)
    // -----------------------------------------------------------------------------------------------

    /**
     * Search the best allele for a read.
     *
     * @param sampleIndex including sample index.
     * @param readIndex  target read index.
     *
     * @return never {@code null}, but with {@link BestAllele#allele allele} == {@code null} if non-could be found.
     */
    private BestAllele searchBestAllele(final int sampleIndex,
                                        final int readIndex,
                                        final boolean canBeReference) {

        final int alleleCount = alleles.numberOfAlleles();

        // if no alleles or
        //    only 1 allele and reference alleles are stored as the index 0 but caller doesn't want ref allele
        // return a un-informative BestAllele with null allele
        if (alleleCount == 0 || (alleleCount == 1 && referenceAlleleIndex == 0 && !canBeReference)) {
            return new BestAllele(sampleIndex, readIndex, -1, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        }

        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        int bestAlleleIndex = canBeReference || referenceAlleleIndex != 0 ? 0 : 1;

        double bestLikelihood = sampleValues[bestAlleleIndex][readIndex];
        double secondBestLikelihood = Double.NEGATIVE_INFINITY;
        for (int a = bestAlleleIndex + 1; a < alleleCount; a++) {
            if (!canBeReference && referenceAlleleIndex == a) {
                continue;
            }
            final double candidateLikelihood = sampleValues[a][readIndex];
            if (candidateLikelihood > bestLikelihood) {
                bestAlleleIndex = a;
                secondBestLikelihood = bestLikelihood;
                bestLikelihood = candidateLikelihood;
            } else if (candidateLikelihood > secondBestLikelihood) {
                secondBestLikelihood = candidateLikelihood;
            }
        }
        return new BestAllele(sampleIndex,readIndex,bestAlleleIndex,bestLikelihood,secondBestLikelihood);
    }

    /**
     * Returns reads stratified by their best allele.
     * @return never {@code null}, perhaps empty.
     */
    @VisibleForTesting
    Map<A,List<GATKRead>> readsByBestAlleleMap() {
        final int alleleCount = alleles.numberOfAlleles();
        final Map<A,List<GATKRead>> result = new LinkedHashMap<>(alleleCount);
        final int totalReadCount = readCount();
        for (int a = 0; a < alleleCount; a++) {
            result.put(alleles.getAllele(a), new ArrayList<>(totalReadCount));
        }
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            readsByBestAlleleMap(s, result);
        }
        return result;
    }

    /**
     * Returns reads stratified by their best allele.
     * @param sampleIndex the target sample.
     * @return never {@code null}, perhaps empty.
     */
    private Map<A,List<GATKRead>> readsByBestAlleleMap(final int sampleIndex) {
        Utils.validIndex(sampleIndex, numberOfSamples());
        final int alleleCount = alleles.numberOfAlleles();
        final int sampleReadCount = readsBySampleIndex[sampleIndex].length;
        final Map<A,List<GATKRead>> result = new LinkedHashMap<>(alleleCount);
        for (int a = 0; a < alleleCount; a++) {
            result.put(alleles.getAllele(a), new ArrayList<>(sampleReadCount));
        }
        readsByBestAlleleMap(sampleIndex,result);
        return result;
    }

    private void readsByBestAlleleMap(final int sampleIndex,
                                      final Map<A,List<GATKRead>> result) {
        final GATKRead[] reads = readsBySampleIndex[sampleIndex];
        final int readCount = reads.length;

        for (int r = 0; r < readCount; r++) {
            final BestAllele bestAllele = searchBestAllele(sampleIndex,r,true);
            if (!bestAllele.isInformative()) {
                continue;
            }
            result.get(bestAllele.allele).add(bestAllele.read);
        }
    }

    /**
     * Returns the index of a read within a sample read-likelihood sub collection.
     * @param sampleIndex   the sample index.
     * @param read          the query read.
     * @return              -1 if there is no such read in that sample, 0 or greater otherwise.
     */
    @VisibleForTesting
    int readIndex(final int sampleIndex,
                  final GATKRead read) {
        final Object2IntMap<GATKRead> readIndex = readIndicesBySampleIndex(sampleIndex);
        return readIndex.containsKey(read) ? readIndex.getInt(read) : -1;
    }

    private Object2IntMap<GATKRead> readIndicesBySampleIndex(final int sampleIndex) {
        if (readIndexBySampleIndex[sampleIndex] == null) {
            final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
            final int sampleReadCount = sampleReads.length;
            readIndexBySampleIndex[sampleIndex] = new Object2IntOpenHashMap<>(sampleReadCount);
            for (int r = 0; r < sampleReadCount; r++) {
                readIndexBySampleIndex[sampleIndex].put(sampleReads[r], r);
            }
        }
        return readIndexBySampleIndex[sampleIndex];
    }

    // -----------------------------------------------------------------------------------------------
    // Trivial accessors
    // -----------------------------------------------------------------------------------------------

    /**
     * Returns the collection of best allele estimates for the reads based on the read-likelihoods.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per read in the read-likelihoods collection.
     */
    public Collection<BestAllele> bestAlleles() {
        final List<BestAllele> result = new ArrayList<>(100); // blind estimate.
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final int readCount = sampleReads.length;
            for (int r = 0; r < readCount; r++) {
                result.add(searchBestAllele(s, r, true));
            }
        }
        return result;
    }

    /**
     * Create an independent copy of this read-likelihoods collection
     */
    public ReadLikelihoods<A> copy() {

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        final double[][][] newLikelihoodValues = new double[sampleCount][alleleCount][];

        @SuppressWarnings({"unchecked", "rawtypes"})
        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];

        for (int s = 0; s < sampleCount; s++) {
            newReadsBySampleIndex[s] = readsBySampleIndex[s].clone();
            for (int a = 0; a < alleleCount; a++) {
                newLikelihoodValues[s][a] = valuesBySampleIndex[s][a].clone();
            }
        }

        // Finally we create the new read-likelihood
        return new ReadLikelihoods<>(
                alleles,
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues);
    }

    /**
     * Returns the index of a sample within the likelihood collection.
     *
     * @param sample the query sample.
     *
     * @throws IllegalArgumentException if {@code sample} is {@code null}.
     * @return -1 if the allele is not included, 0 or greater otherwise.
     */
    public int indexOfSample(final String sample) {
        return samples.indexOfSample(sample);
    }

    /**
     * Number of samples included in the likelihood collection.
     * @return 0 or greater.
     */
    public int numberOfSamples() {
        return samples.numberOfSamples();
    }

    /**
     * Returns sample name given its index.
     *
     * @param sampleIndex query index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is negative.
     *
     * @return never {@code null}.
     */
    public String getSample(final int sampleIndex) {
        return samples.getSample(sampleIndex);
    }

    /**
     * Returns the index of an allele within the likelihood collection.
     *
     * @param allele the query allele.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     *
     * @return -1 if the allele is not included, 0 or greater otherwise.
     */
    public int indexOfAllele(final A allele) {
        return alleles.indexOfAllele(allele);
    }

    /**
     * Returns number of alleles in the collection.
     * @return 0 or greater.
     */
    public int numberOfAlleles() {
        return alleles.numberOfAlleles();
    }

    /**
     * Returns the allele given its index.
     *
     * @param alleleIndex the allele index.
     *
     * @throws IllegalArgumentException the allele index is {@code null}.
     *
     * @return never {@code null}.
     */
    public A getAllele(final int alleleIndex) {
        return alleles.getAllele(alleleIndex);
    }

    /**
     * Returns the samples in this read-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable view on the read-likelihoods sample list.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<String> samples() {
        return Collections.unmodifiableList(sampleList == null ? sampleList = samples.asListOfSamples() : sampleList);
    }

    /**
     * HYQ_doc_log: documentation is an apparent copy-paste error.
     * Returns the samples in this read-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable. It will not be updated if the collection
     *     allele list changes.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<A> alleles() {
        return Collections.unmodifiableList(alleleList == null ? alleleList = alleles.asListOfAlleles() : alleleList);
    }

    /**
     * Returns the reads that belong to a sample sorted by their index (within that sample).
     *
     * @param sampleIndex   the requested sample.
     * @return              never {@code null} but perhaps a zero-length array if there is no reads in sample.
     *                      No element in the array will be {@code null}.
     */
    public List<GATKRead> sampleReads(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        final List<GATKRead> extantList = readListBySampleIndex[sampleIndex];
        if (extantList == null) {
            return readListBySampleIndex[sampleIndex] = Collections.unmodifiableList(Arrays.asList(readsBySampleIndex[sampleIndex]));
        } else {
            return extantList;
        }
    }

    /**
     * Returns a read vs allele likelihood matrix corresponding to a sample.
     *
     * @param sampleIndex target sample.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not null.
     *
     * @return never {@code null}
     */
    public LikelihoodMatrix<A> sampleMatrix(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        final LikelihoodMatrix<A> extantResult = sampleMatrices[sampleIndex];
        if (extantResult == null) {
            return sampleMatrices[sampleIndex] = new SampleMatrix(sampleIndex);
        } else {
            return extantResult;
        }
    }

    /**
     * Returns the total number of reads in the read-likelihood collection.
     *
     * @return never {@code null}
     */
    public int readCount() {
        int sum = 0;
        final int sampleCount = samples.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            sum += readsBySampleIndex[i].length;
        }
        return sum;
    }

    /**
     * Returns the number of reads that belong to a sample in the read-likelihood collection.
     * @param sampleIndex the query sample index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not a valid sample index.
     * @return 0 or greater.
     */
    public int sampleReadCount(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        return readsBySampleIndex[sampleIndex].length;
    }

    // -----------------------------------------------------------------------------------------------
    // Like a tool
    // TODO: these should be refactored to utilities. they don't depend on data stored here.
    // -----------------------------------------------------------------------------------------------

    public static boolean unclippedReadOverlapsRegion(final GATKRead read,
                                                      final Locatable region) {

        final String readReference = read.getContig();
        if (!Objects.equals(readReference, region.getContig())) {
            return false;
        }

        final int readStart = read.getUnclippedStart();
        if (readStart > region.getEnd()) {
            return false;
        }

        final int readEnd = read.isUnmapped() ? read.getUnclippedEnd() : Math.max(read.getUnclippedEnd(), read.getUnclippedStart());
        return readEnd >= region.getStart();
    }

    /**
     * Given a collection of likelihood in the old map format, it creates the corresponding read-likelihoods collection.
     *
     * @param map the likelihoods to transform.
     *
     * @throws IllegalArgumentException if {@code map} is {@code null}.
     *
     * @return never {@code null}.
     */
    public static ReadLikelihoods<Allele> fromPerAlleleReadLikelihoodsMap(final Map<String,PerReadAlleleLikelihoodMap> map) {
        Utils.nonNull(map);

        // First we need to create the read-likelihood collection with all required alleles, samples and reads.
        final SampleList sampleList = new IndexedSampleList(map.keySet());
        final Set<Allele> alleles = new LinkedHashSet<>(10);
        final Map<String,List<GATKRead>> sampleToReads = new LinkedHashMap<>(sampleList.numberOfSamples());
        for (final Map.Entry<String,PerReadAlleleLikelihoodMap> entry : map.entrySet()) {
            final String sample = entry.getKey();
            final PerReadAlleleLikelihoodMap sampleLikelihoods = entry.getValue();
            alleles.addAll(sampleLikelihoods.getAllelesSet());
            sampleToReads.put(sample,new ArrayList<>(sampleLikelihoods.getLikelihoodReadMap().keySet()));
        }

        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
        final ReadLikelihoods<Allele> result = new ReadLikelihoods<>(sampleList,alleleList,sampleToReads);

        // Now set the likelihoods.
        for (final Map.Entry<String,PerReadAlleleLikelihoodMap> sampleEntry : map.entrySet()) {
            final LikelihoodMatrix<Allele> sampleMatrix = result.sampleMatrix(result.indexOfSample(sampleEntry.getKey()));
            for (final Map.Entry<GATKRead,Map<Allele,Double>> readEntry : sampleEntry.getValue().getLikelihoodReadMap().entrySet()) {
                final GATKRead read = readEntry.getKey();
                final int readIndex = sampleMatrix.indexOfRead(read);
                for (final Map.Entry<Allele,Double> alleleEntry : readEntry.getValue().entrySet()) {
                    final int alleleIndex = result.indexOfAllele(alleleEntry.getKey());
                    sampleMatrix.set(alleleIndex,readIndex,alleleEntry.getValue());
                }
            }
        }
        return result;
    }

    /**
     * Transform into a multi-sample HashMap backed {@link PerReadAlleleLikelihoodMap} type.
     * @return never {@code null}.
     *
     */
    public Map<String, PerReadAlleleLikelihoodMap> toPerReadAlleleLikelihoodMap() {
        final int sampleCount = samples.numberOfSamples();
        final Map<String, PerReadAlleleLikelihoodMap> result = new LinkedHashMap<>(sampleCount);
        for (int s = 0; s < sampleCount; s++)
            result.put(samples.getSample(s),toPerReadAlleleLikelihoodMap(s));
        return result;
    }

    /**
     * Transform into a single-sample HashMap backed {@link PerReadAlleleLikelihoodMap} type.
     *
     * @return never {@code null}.
     */
    public PerReadAlleleLikelihoodMap toPerReadAlleleLikelihoodMap(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        final PerReadAlleleLikelihoodMap result = new PerReadAlleleLikelihoodMap();
        final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];

        final int alleleCount = alleles.numberOfAlleles();
        final int sampleReadCount = sampleReads.length;
        for (int a = 0; a < alleleCount; a++) {
            final A allele = alleles.getAllele(a);
            final double[] readLikelihoods = valuesBySampleIndex[sampleIndex][a];
            for (int r = 0; r < sampleReadCount; r++)
                result.add(sampleReads[r], allele, readLikelihoods[r]);
        }
        return result;
    }
}