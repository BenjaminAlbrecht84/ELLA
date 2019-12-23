package ella.model.aligner.aligning.extending;

import java.util.ArrayList;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.aligning.query.ReferenceHit;
import ella.model.aligner.aligning.seeding.Seed;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.CodonTranslator;
import ella.model.aligner.utils.AlignmentCompressor;
import ella.model.aligner.utils.MyStrand;
import ella.model.aligner.utils.frameshiftAligner.Aligner_Mod9;
import ella.model.aligner.utils.frameshiftAligner.DP_Matrix;
import ella.model.aligner.utils.frameshiftAligner.DP_Matrix_last;
import ella.model.aligner.utils.frameshiftAligner.Zheng_XDrop_Frameshift.AliMode;
import ella.model.aligner.utils.gotoAligner.GotohAligner;
import ella.model.aligner.utils.gotoAligner.GotohAligner_Mod1;
import ella.model.aligner.utils.gotoAligner.GotohAligner_Mod2;
import ella.model.io.MyParameters;

public class AlignmentExtender {

    //	private Zheng_XDrop_Frameshift frameShiftAligner = new Zheng_XDrop_Frameshift();
//	private Aligner_Mod1 modOneAligner = new Aligner_Mod1();
//	private Aligner_Mod2 modTwoAligner = new Aligner_Mod2();
//	private Aligner_Mod3 modThreeAligner = new Aligner_Mod3();
//	private Aligner_Mod4 modFourAligner = new Aligner_Mod4();
//	private Aligner_Mod5 modFiveAligner = new Aligner_Mod5();
//	private Aligner_Mod6 modSixAligner = new Aligner_Mod6();
//	private Aligner_Mod7 modSevenAligner = new Aligner_Mod7();
//	private Aligner_Mod8 modEightAligner = new Aligner_Mod8();
    //	private Aligner_Mod10 modTenAligner = new Aligner_Mod10();
//	private Aligner_Mod11 modElevenAligner = new Aligner_Mod11();
//	private Aligner_Mod13 modThirteenAligner = new Aligner_Mod13();

    private GotohAligner gotohAligner = new GotohAligner();
    private GotohAligner_Mod1 gotohAligner1 = new GotohAligner_Mod1();
    private final GotohAligner_Mod2 gotohAligner2 = new GotohAligner_Mod2();
    private final Aligner_Mod9 modNineAligner = new Aligner_Mod9();
    private Aligner aligner;

    public void run(ReferenceHit refHit, IndexText text, String query, QueryContainer qC, DP_Matrix dpMatrix, Aligner aligner,
                    ArrayList<FrameRun> resultingRuns, boolean b) {

        this.aligner = aligner;
        ArrayList<Seed> seeds = refHit.getSeedChain().getSeeds();
        ArrayList<FrameHit> frameHits = new ArrayList<>();
        int firstSeedFrame = seeds.get(0).getFrame().getNumericalID();

        // extending first seed to the left end
        Seed firstSeed = seeds.get(0);
        int queryLength = firstSeed.getQueryStart();
        int refLength = firstSeed.getRefStart();
        int expQueryLength = (int) Math.round(new Double(refLength * 3) * 0.97);
        int expRefLength = (int) Math.round(new Double(queryLength / 3) * 1.03);
        int queryStartOffset = 0, refStartOffset = 0;
        if (queryLength > expQueryLength && MyParameters.FRAMESHIFT_PENALTY != null)
            queryStartOffset = queryLength - expQueryLength;
        else if (refLength > expRefLength && MyParameters.FRAMESHIFT_PENALTY != null)
            refStartOffset = refLength - expRefLength;
        int refStart = refStartOffset;
        int queryStart = 0 + queryStartOffset;
        int refEnd = firstSeed.getRefStart();
        int queryEnd = firstSeed.getQueryStart() + 1;
        Object[] aliResult = cmpAlignment(refHit.getTextReferenceStart(), refStart, refEnd, queryStart, queryEnd, text, query, AliMode.LEFT,
                dpMatrix, firstSeedFrame);
        frameHits.addAll(
                generateFrameHits(aliResult, queryStart + (int) aliResult[4] * 3, refStart + (int) aliResult[5], refHit.getStrand(), query.length()));
        frameHits.addAll(generateFrameHits(frameHits.size() != 0 ? frameHits.get(frameHits.size() - 1) : null, seeds.get(0), text, query));

        // extending seeds to the respective seed on the right
        int i;
        for (i = 0; i < seeds.size() - 1; i++) {

            Seed leftSeed = seeds.get(i);
            Seed rightSeed = seeds.get(i + 1);
            refStart = leftSeed.getRefStart() + leftSeed.getRefLength();
            queryStart = leftSeed.getQueryStart() + leftSeed.getQueryLength() - 1;
            refEnd = rightSeed.getRefStart();
            queryEnd = rightSeed.getQueryStart() + 1;
            queryStart = MyParameters.FRAMESHIFT_PENALTY != null ? queryStart : queryStart + 1;
            aliResult = cmpAlignment(refHit.getTextReferenceStart(), refStart, refEnd, queryStart, queryEnd, text, query, AliMode.MIDDLE, dpMatrix, firstSeedFrame);

            // checking if xDrop-Alignment closed the gap to the right neighbored seed
            if ((int) aliResult[5] != 0) {
                ReferenceHit subRefHit = refHit.extractSubHit(i + 1);
                run(subRefHit, text, query, qC, dpMatrix, aligner, resultingRuns, true);
                break;
            }

            frameHits.addAll(generateFrameHits(aliResult, queryStart, refStart, refHit.getStrand(), query.length()));
            frameHits.addAll(generateFrameHits(frameHits.get(frameHits.size() - 1), seeds.get(i + 1), text, query));

        }

        // extending last seed to the right end
        Seed lastSeed = seeds.get(i);
        queryLength = query.length() - lastSeed.getQueryEnd() + 1;
        refLength = (int) (refHit.getTextReferenceStart() + refHit.getTextReferenceLength() - lastSeed.getIndexEnd());
        expQueryLength = (int) Math.round(new Double(refLength * 3) * 0.97);
        expRefLength = (int) Math.round(new Double(queryLength / 3) * 1.03);
        refStart = lastSeed.getRefStart() + lastSeed.getRefLength();
        queryStart = lastSeed.getQueryStart() + lastSeed.getQueryLength() - 1;
        refEnd = refStart + Math.min(refLength, expRefLength);
        queryEnd = queryStart + Math.min(queryLength, expQueryLength);
        queryStart = MyParameters.FRAMESHIFT_PENALTY != null ? queryStart : queryStart + 1;
        aliResult = cmpAlignment(refHit.getTextReferenceStart(), refStart, refEnd, queryStart, queryEnd, text, query, AliMode.RIGHT, dpMatrix, firstSeedFrame);
        frameHits.addAll(generateFrameHits(aliResult, queryStart, refStart, refHit.getStrand(), query.length()));

        // merging neighboring frame-hits of the same frame
        frameHits = mergeFrameHits(frameHits);
        FrameRun frameRun = new FrameRun(frameHits, refHit, qC, text);

        // checking if alignment is optimal
        if (!frameRun.getFrameShiftHit().isOptimal(query))
            return;

        // // checking eValue, coverage, and raw score
        if (frameRun.getEValue() > MyParameters.MAX_EVALUE || frameRun.getRefCoverage() < MyParameters.MIN_COVERAGE || frameRun.getRawScore() < MyParameters.X_BESTSCORE_DROP)
            return;

        resultingRuns.add(frameRun);

    }

    private ArrayList<FrameHit> mergeFrameHits(ArrayList<FrameHit> frameHits) {
        ArrayList<FrameHit> mergedFrameHits = new ArrayList<FrameHit>();
        mergedFrameHits.add(frameHits.get(0));
        for (int i = 1; i < frameHits.size(); i++) {
            FrameHit f1 = mergedFrameHits.get(mergedFrameHits.size() - 1);
            FrameHit f2 = frameHits.get(i);
            if (f1.getFrame() == f2.getFrame() && f1.getRefEnd() == f2.getRefStart() && f1.getQueryEnd() == f2.getQueryStart()) {
                FrameHit f = f1.addFrameHit(f2);
                mergedFrameHits.remove(mergedFrameHits.size() - 1);
                mergedFrameHits.add(f);
            } else
                mergedFrameHits.add(f2);
        }
        return mergedFrameHits;
    }

    private Object[] cmpAlignment(long textRefStart, int refStart, int refEnd, int queryStart, int queryEnd, IndexText text, String query,
                                  AliMode aliMode, DP_Matrix dpMatrix, int frame) {

        // System.out.println("\n[" + refStart + "," + refEnd + "] ");
        // System.out.println("[" + queryStart + "," + queryEnd + "] ");

        long[] proteinCoord = {textRefStart + refStart, textRefStart + refEnd};
        String refSeq = extractRefSeq(textRefStart + refStart, textRefStart + refEnd, text);
        String querySeq = query.substring(queryStart, queryEnd);

        // System.out.println(querySeq);
        // System.out.println(refSeq);

        Object[] result;
        if (MyParameters.FRAMESHIFT_PENALTY != null && MyParameters.FRAMESHIFT_PENALTY > 0) {
            // Object[] result = frameShiftAligner.run(querySeq, refSeq, aliMode, MyParameters.X_DIAGONAL_DROP);
            // Object[] result = modOneAligner.run(querySeq, refSeq, aliMode, MyParameters.X_DIAGONAL_DROP, null);
            // Object[] result = modTwoAligner.run(querySeq, refSeq, aliMode, MyParameters.X_DIAGONAL_DROP, null);
            // Object[] result = modThreeAligner.run(querySeq, refSeq, aliMode, MyParameters.X_DIAGONAL_DROP, null);
            // Object[] result = modFourAligner.run(querySeq, refSeq, aliMode, MyParameters.X_DIAGONAL_DROP, null);
            // Object[] result = modFiveAligner.run(querySeq, refSeq, aliMode, MyParameters.X_DIAGONAL_DROP, dpMatrix);
            // Object[] result = modSixAligner.run(querySeq, refSeq, aliMode, MyParameters.X_DIAGONAL_DROP, dpMatrix);
            // Object[] result = modSevenAligner.run(querySeq, refSeq, aliMode, MyParameters.X_DIAGONAL_DROP, null);
            // Object[] result = modEightAligner.run(querySeq, refSeq, aliMode, null);
            result = modNineAligner.run(querySeq, refSeq, aliMode, dpMatrix, aligner);
            // Object[] result = modTenAligner.run(querySeq, aliMode, null, proteinCoord, text);
            // Object[] result = modElevenAligner.run(querySeq, refSeq, aliMode, null, ella.model.aligner);
            // Object[] result = modThirteenAligner.run(querySeq, refSeq, aliMode, dpMatrix);
        } else {
//            result = gotohAligner.run(querySeq, refSeq, frame, aliMode);
//            result = gotohAligner1.run(querySeq, refSeq, frame, aliMode);
            result = gotohAligner2.run(querySeq, refSeq, frame, aliMode);
        }

        // System.out.println("\n" + result[0] + "\n" + result[1] + "\n" + result[2]);

        return result;
    }

    public String extractRefSeq(long indexTextStart, long indexTextEnd, IndexText text) {

        if (indexTextEnd - indexTextStart < 0)
            System.exit(0);

        StringBuffer buf = new StringBuffer((int) (indexTextEnd - indexTextStart + 1));
        long i = indexTextStart;
        int c;
        while (i < indexTextEnd && (c = text.readPosition(i)) != 127) {
            buf.append(Alphabet.getCharacter(c));
            i++;
        }
        return buf.toString();
    }

    private ArrayList<FrameHit> generateFrameHits(FrameHit lastAddedHit, Seed seed, IndexText text, String query) {
        double positionDiff = lastAddedHit != null ? seed.getQueryStart() - lastAddedHit.getQueryEnd() : 0;
        int offset = positionDiff > 1 ? (int) Math.ceil((positionDiff - 1) / 3.) * 3 : 0;
        offset = positionDiff < -1 ? (int) Math.ceil((positionDiff + 1) / 3.) * 3 : offset;

        String queryAli = CodonTranslator.translate(query.substring(seed.getQueryStart() - offset, seed.getQueryEnd()), 1);
        String qGapString = offset < 0 ? generateGapString(-offset / 3) : "";
        queryAli = qGapString + queryAli;
        String refAli = extractRefSeq(seed.getIndexStart(), seed.getIndexEnd(), text);
        refAli = offset > 0 ? generateGapString(offset / 3) + refAli : refAli;

        String frameAli = generateFrameString(seed, refAli.length());
        String[] aliResult = {queryAli, refAli, frameAli};
        ArrayList<FrameHit> frameHits = new ArrayList<FrameHit>();

        frameHits.add(setUpFrameHit(aliResult, seed.getQueryStart() - offset, seed.getRefStart(), (queryAli.length() - qGapString.length()) * 3,
                seed.getRefLength(), seed.getStrand(), query.length()));

        // System.out.println("\n" + queryAli + "\n" + refAli + "\n");

        return frameHits;
    }

    private String generateFrameString(Seed seed, int l) {
        int frame = Math.abs(seed.getFrame().getNumericalID());
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < l; i++)
            buf.append(frame);
        return buf.toString();
    }

    private String generateGapString(int l) {
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < l; i++)
            buf.append("-");
        return buf.toString();
    }

    private ArrayList<FrameHit> generateFrameHits(Object[] aliResult, int queryStart, int refStart, MyStrand strand, int totalQueryLength) {

        ArrayList<FrameHit> frameHits = new ArrayList<FrameHit>();

        StringBuilder[] frameAliResult = {new StringBuilder(), new StringBuilder(), new StringBuilder()};
        int qStart = queryStart, qLength = 0, rStart = refStart, rLength = 0;

        String frameIDs = (String) aliResult[2];
        if (frameIDs.isEmpty())
            return frameHits;

        Integer lastFrameID = 1;
        for (int i = 0; i <= frameIDs.length(); i++) {

            Integer frameID = i < frameIDs.length() ? Character.getNumericValue(frameIDs.charAt(i)) : null;

            if (i == 0)
                qStart = qStart + frameID - 1;

            if (i < frameIDs.length() && (i == 0 || lastFrameID == frameID)) { // extending frame hit
                char q = ((String) aliResult[0]).charAt(i);
                char r = ((String) aliResult[1]).charAt(i);
                frameAliResult[0] = frameAliResult[0].append(q);
                frameAliResult[1] = frameAliResult[1].append(r);
                qLength = q != '-' ? qLength + 3 : qLength;
                rLength = r != '-' ? rLength + 1 : rLength;
            } else { // reporting frame hit

                String queryAli = frameAliResult[0].toString();
                String refAli = frameAliResult[1].toString();
                String[] subAli = {queryAli, refAli};
                if (rLength >= 1) {
                    FrameHit h = setUpFrameHit(subAli, qStart, rStart, qLength, rLength, strand, totalQueryLength);
                    frameHits.add(h);
                }

                if (i < frameIDs.length()) {

                    // computing offset at query start
                    int offset = frameID - lastFrameID;

                    // artifact produced by the Zheng-Frameshift-Aligner
                    if (i > 0) {
                        switch (offset) {
                            case (-2):
                                offset = 1;
                                break;
                            case (-1):
                                offset = -1;
                                break;
                            case (2):
                                offset = -1;
                                break;
                            case (1):
                                offset = 1;
                                break;
                        }
                    }

                    // resetting parameters for recording next frame hit
                    qStart = qStart + qLength + offset;
                    rStart = rStart + rLength;
                    qLength = 0;
                    rLength = 0;
                    frameAliResult[0] = new StringBuilder();
                    frameAliResult[1] = new StringBuilder();
                    i--;

                }

            }

            lastFrameID = frameID;

        }

        return frameHits;

    }

    private FrameHit setUpFrameHit(String[] aliResult, int queryStart, int refStart, int qLength, int rLength, MyStrand strand,
                                   int totalQueryLength) {
        int rawScore = MyParameters.SCORING_MATRIX.cmpAlignmentScore(aliResult[0], aliResult[1]);
        int frame = strand == MyStrand.FWD ? (queryStart % 3) + 1 : -((queryStart % 3) + 1);
        ArrayList<Byte> editOperations = AlignmentCompressor.run(aliResult);
        return new FrameHit(frame, rawScore, queryStart, refStart, qLength, rLength, editOperations);
    }

}
