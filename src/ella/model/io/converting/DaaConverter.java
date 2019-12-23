package ella.model.io.converting;

import ella.model.aligner.utils.BlastStatisticsHelper;
import ella.model.aligner.utils.ScoringMatrix;
import ella.presenter.Presenter;

import javax.xml.transform.TransformerException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class DaaConverter {

    public static String[] BLAST_FORMATS = {"BLAST Tab (*.tab)", "BLAST pairwise (*.txt)"};

    public long n;
    public double K, LAMBDA, lnK, LN2 = (float) Math.log(2);
    private DaaReader daaReader;
    private Presenter presenter;

    public void run(Presenter presenter, String daaFilePath, String outFilePath, String format, String ncpus) {
        this.presenter = presenter;
        File daaFile = new File(daaFilePath);
        Integer cores = Integer.parseInt(ncpus);
        if (format.equals(BLAST_FORMATS[0]))
            toBlastTabFormat(daaFile, cores, outFilePath, true);
        else if (format.equals(BLAST_FORMATS[1]))
            toBlastPairwise(daaFile, cores, outFilePath);
        presenter.reportConvertStatus("File '" + daaFile.getName() + "' successfully converted!", false);
    }

    public void toBlastPairwise(File daaFile, int cores, String outFile) {
        daaReader = new DaaReader(daaFile, false);
        daaReader.totalProgressProperty().addListener((a, b, c) -> reportProgress("Reading DAA File... (" + c + "%)"));
        initializeStatistics();
        daaReader.parseAllHits(cores);
        new BlastPairwiseFile().write(daaReader, outFile, this);
    }

    public void toBlastXMLFormat(File daaFile, int cores, String outFile) {

        daaReader = new DaaReader(daaFile, false);
        daaReader.totalProgressProperty().addListener((a, b, c) -> reportProgress("Reading DAA File... (" + c + "%)"));
        initializeStatistics();
        ArrayList<DaaHit> allHits = daaReader.parseAllHits(cores);

        BlastXMLFile blastXMLFile = new BlastXMLFile();
        blastXMLFile.initializeFile("blastx", "ELLA - Enhanced Local Aligner", "UNKNOWN", "UNKNOWN", "Query_1",
                allHits.get(0).getQueryName(), String.valueOf(allHits.get(0).getTotalQueryLength()));
        blastXMLFile.addUsedProgramParameters(daaReader.getDaaHeader().getScoreMatrix(), "0",
                String.valueOf(daaReader.getDaaHeader().getGapOpen()), String.valueOf(daaReader.getDaaHeader().getGapExtend()),
                "F");

        int iteration = 1;
        for (String queryId : daaReader.getReadId2Hits().keySet()) {
            blastXMLFile.addHitIteration(iteration, daaReader.getReadId2Hits().get(queryId), daaReader.getDaaHeader().getDbLetters(), daaReader.getDaaHeader().getDbSeqs(), daaReader, this);
            iteration++;
        }

        try {
            blastXMLFile.writeXML(outFile);
        } catch (TransformerException e) {
            e.printStackTrace();
        }
    }

    public void toBlastTabFormat(File daaFile, int cores, String outFile, boolean writeHeaderLine) {
        daaReader = new DaaReader(daaFile, false);
        daaReader.totalProgressProperty().addListener((a, b, c) -> reportProgress("Reading DAA File... (" + c + "%)"));
        initializeStatistics();
        new BlastTabFile().write(daaReader, outFile, cores, writeHeaderLine, this);
    }

    public void initializeStatistics() {
        String matrix = daaReader.getDaaHeader().getScoreMatrix();
        int gapOpen = daaReader.getDaaHeader().getGapOpen();
        int gapExtend = daaReader.getDaaHeader().getGapExtend();
        long dbLetters = daaReader.getDaaHeader().getDbLetters().longValue();
        ScoringMatrix scoringMatrix = new ScoringMatrix(matrix, gapOpen, gapExtend);
        double[] stats = BlastStatisticsHelper.getStatistics(scoringMatrix);
        LAMBDA = stats[0];
        K = stats[1];
        n = dbLetters;
        lnK = Math.log(K);
    }

    public double getEValue(int alignmentScore, int queryLength) {
        return K * n * (double) queryLength * Math.exp(-LAMBDA * alignmentScore);
    }

    public double getBitScore(int alignmentScore) {
        return (LAMBDA * alignmentScore - lnK) / LN2;
    }

    public void reportProgress(String message) {
        if (presenter != null && daaReader != null) {
            String daaFile = daaReader.getDaaFile().getName();
            presenter.reportConvertStatus(daaFile + ": " + message, true);
        }
    }
}
