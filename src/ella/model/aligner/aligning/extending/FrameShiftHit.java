package ella.model.aligner.aligning.extending;

import java.util.ArrayList;

import ella.model.aligner.utils.AlignmentCompressor;
import ella.model.aligner.utils.AlignmentDecompressor;
import ella.model.aligner.utils.BlastStatisticsHelper;
import ella.model.io.MyParameters;

public class FrameShiftHit {

    private int frame;
    private int rawScore;
    private int queryStart, refStart;
    private ArrayList<Byte> editOperations;

    public FrameShiftHit(ArrayList<FrameHit> hits, long textReferenceStart) {

        FrameHit h0 = hits.get(0);
        this.frame = h0.getFrame();
        this.queryStart = h0.getQueryStart();
        this.refStart = h0.getRefStart();

        // computing rawScore and frameshift alignment operations
        editOperations = h0.getEditOperations();
        rawScore = h0.getScore();
        for (int i = 1; i < hits.size(); i++) {
            int queryDiff = hits.get(i).getQueryStart() - hits.get(i - 1).getQueryEnd();
            if (MyParameters.FRAMESHIFT_PENALTY != null)
                rawScore -= Math.abs(queryDiff) * MyParameters.FRAMESHIFT_PENALTY;
            rawScore += hits.get(i).getScore();
            editOperations.addAll(AlignmentCompressor.cmpFrameShiftOperations(queryDiff));
            editOperations.addAll(hits.get(i).getEditOperations());
        }

    }

    public boolean isOptimal(String queryDNA) {
        return AlignmentDecompressor.isOptimalAlignment(queryDNA, queryStart, editOperations);
    }

    public int getRawScore() {
        return rawScore;
    }

    public double getBitScore() {
        return BlastStatisticsHelper.getBitScore(rawScore);
    }

    public int getFrame() {
        return frame;
    }

    public int getQueryStart() {
        return queryStart;
    }

    public int getRefStart() {
        return refStart;
    }

    public ArrayList<Byte> getEditOperations() {
        return editOperations;
    }

    public double cmpEValue(int totalQueryLength) {
        return BlastStatisticsHelper.getEValue(rawScore, totalQueryLength);
    }

    public double cmpEG2Value() {
        return BlastStatisticsHelper.getEG2Value(rawScore);
    }

    public String[] getAlignment(String queryDNA) {
        return AlignmentDecompressor.getAlignment(queryDNA, queryStart, editOperations);
    }

    public void printAlignment(String query) {
        String ali[] = AlignmentDecompressor.getAlignment(query, queryStart, editOperations);
        System.out.println(ali[0] + "\n" + ali[1] + "\n" + ali[2]);
    }

}
