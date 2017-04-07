package org.broadinstitute.hellbender.tools;

import com.google.api.services.genomics.Genomics;
import com.intel.genomicsdb.GenomicsDBExportConfiguration;
import com.intel.genomicsdb.GenomicsDBFeatureReader;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.io.FileUtils;
>>>>>>> 25e176b... genomicsdb importer test is ready
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public final class GenomicsDBImportIntegrationTest extends CommandLineProgramTest {

  private static final File TEST_DATA_DIR = getTestDataDir();
  private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/genomicsdb";
  private static final File GENOMICSDB_WORKSPACE = new File(TEST_OUTPUT_DIRECTORY + "/tiledb-ws");
  private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/genomicsdb";
  private static final File GENOMICSDB_WORKSPACE = new File(TEST_OUTPUT_DIRECTORY + "/tiledb-ws");
  private static final String GENOMICSDB_ARRAYNAME = "gatk4-genomicsdb-test-0";
  private static final File TEST_CALLSETMAP_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/callset.json");
  private static final File TEST_VIDMAP_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/vidmap.json");
  private static final File TEST_LOADER_JSON_FILE = new File(TEST_OUTPUT_DIRECTORY + "/loader.json");
  private static final String TEST_GENERATED_COMBINED_GVCF = TEST_OUTPUT_DIRECTORY + "/test.combined.g.vcf.gz";

  private static final String hg00096 = publicTestDir + "large/gvcfs/HG00096.g.vcf.gz";
  private static final String hg00268 = publicTestDir + "large/gvcfs/HG00268.g.vcf.gz";
  private static final String na19625 = publicTestDir + "large/gvcfs/NA19625.g.vcf.gz";
  private static final String combined = publicTestDir + "large/gvcfs/combined.gatk3.7.g.vcf.gz";

  private static final File TEST_REFERENCE_GENOME = new File(publicTestDir + "/large/Homo_sapiens_assembly38.20.21.fasta");

  @Override
  public String getTestedClassName() {
    return GenomicsDBImport.class.getSimpleName();
  }

  @Test
  public void testGenomicsDBImporter() throws IOException {

    final ArgumentsBuilder args = new ArgumentsBuilder();
    args.add("-GW"); args.add(GENOMICSDB_WORKSPACE.getAbsolutePath());

    SimpleInterval simpleInterval = new SimpleInterval("chr20", 69491, 69521);
    args.add("-L"); args.add(simpleInterval);
    args.add("-V"); args.add(hg00096);
    args.add("-V"); args.add(hg00268);
    args.add("-V"); args.add(na19625);
    args.add("-GA"); args.add(GENOMICSDB_ARRAYNAME);
    args.add("-GVID"); args.add(TEST_VIDMAP_JSON_FILE.getAbsolutePath());
    args.add("-GCS"); args.add(TEST_CALLSETMAP_JSON_FILE.getAbsolutePath());
    args.add("-R"); args.add(TEST_REFERENCE_GENOME);

    runCommandLine(args);

    GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> featureReader =
      new GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream>(
        TEST_VIDMAP_JSON_FILE.getAbsolutePath(),
        TEST_CALLSETMAP_JSON_FILE.getAbsolutePath(),
        GENOMICSDB_WORKSPACE.getAbsolutePath(),
        GENOMICSDB_ARRAYNAME,
        TEST_REFERENCE_GENOME.getAbsolutePath(), null, new BCF2Codec());

    final VariantContextWriter writer =
      new VariantContextWriterBuilder().setOutputVCFStream(new FileOutputStream(TEST_GENERATED_COMBINED_GVCF)).unsetOption(
        Options.INDEX_ON_THE_FLY).build();

    writer.writeHeader((VCFHeader)(featureReader.getHeader()));

    CloseableTribbleIterator<VariantContext> actualVcs = featureReader.query(simpleInterval.getContig(), simpleInterval.getStart(), simpleInterval.getEnd());
    while (actualVcs.hasNext()) {
      writer.add(actualVcs.next());
    }

    actualVcs.close();

    AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(combined, new VCFCodec(), false);
    CloseableTribbleIterator<VariantContext> expectedVcs = reader.query(simpleInterval.getContig(), simpleInterval.getStart(), simpleInterval.getEnd());

    assertCondition(actualVcs, expectedVcs, (a,e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a,e, Collections.emptyList()));

    IOUtils.deleteRecursivelyOnExit(GENOMICSDB_WORKSPACE);
    FileUtils.deleteQuietly(TEST_CALLSETMAP_JSON_FILE);
    FileUtils.deleteQuietly(TEST_VIDMAP_JSON_FILE);
    FileUtils.deleteQuietly(TEST_LOADER_JSON_FILE);
  }

  private static <T> void assertCondition(Iterable<T> actual, Iterable<T> expected, BiConsumer<T,T> assertion){
    final Iterator<T> iterActual = actual.iterator();
    final Iterator<T> iterExpected = expected.iterator();
    while(iterActual.hasNext() && iterExpected.hasNext()){
      assertion.accept(iterActual.next(), iterExpected.next());
    }
    if (iterActual.hasNext()){
      Assert.fail("actual is longer than expected with at least one additional element: " + iterActual.next());
    }
    if (iterExpected.hasNext()){
      Assert.fail("actual is shorter than expected, missing at least one element: " + iterExpected.next());
    }
  }
}
