package org.broadinstitute.hellbender.tools.spark.sv;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 * Represents a collection of relevant information on one locally assembled contig from short reads.
 */
class LocalAssemblyContig implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * This contig was assembled based on FASTQ records collected by {@link FindBadGenomicKmersSpark} and numbered by it.
     */
    final long assemblyID;

    /**
     * This together with {@link #assemblyID} completely determines which contig we are talking about.
     */
    final String contigID;

    final String seq;

    /**
     * Per base coverage for the assembled contig
     * May be {@code null} if assembler doesn't provide this feature (e.g. SGA).
     */
    final byte[] perBaseCoverage;

    /**
     * Read count for the contigs.
     * May be {@code null} if assembler doesn't provide this feature (e.g. SGA).
     */
    final Integer readCountSupport;

    private List<AlignmentRegion> alignmentRecords;

    /**
     * Construct a contig with assembly id, contig/vertex id, bases, but without number of reads generating the contig, or qual.
     */
    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq){
        this(assemblyID, contigID, seq, null, null, null);
    }

    /**
     * Construct a contig with assembly id, contig/vertex id, bases, number of reads generating the contig, and qual.
     * Currently fermi-lite supports the last two features, not SGA.
     */
    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq,
                        final List<AlignmentRegion> alignmentRegions,
                        final Integer readCountSupport, final byte[] perBaseCoverage){
        this.assemblyID = assemblyID;
        this.contigID   = contigID;
        this.seq        = seq;
        this.alignmentRecords = alignmentRegions;
        this.perBaseCoverage = (perBaseCoverage == null) ? null : perBaseCoverage;
        this.readCountSupport = readCountSupport;
    }

    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq,
                        final List<AlignmentRegion> alignmentRegions){
        this(assemblyID, contigID, seq, alignmentRegions, null, null);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        LocalAssemblyContig that = (LocalAssemblyContig) o;

        if (assemblyID != that.assemblyID) return false;
        if (!contigID.equals(that.contigID)) return false;
        if (seq.equals(that.seq)) return false;

        return (perBaseCoverage !=null && that.perBaseCoverage !=null) ? Arrays.equals(perBaseCoverage, that.perBaseCoverage) : (perBaseCoverage ==null && that.perBaseCoverage ==null);
    }

    @Override
    public int hashCode() {
        int result = (int) (assemblyID ^ (assemblyID >>> 32));
        result = 31 * result + contigID.hashCode();
        result = 31 * result + seq.hashCode();
        return result;
    }
}
