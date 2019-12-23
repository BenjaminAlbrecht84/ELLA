package ella.model.io.converting;

import ella.model.aligner.utils.ScoringMatrix;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.ConcurrentHashMap;

public class BlastPairwiseFile {

    private BufferedWriter writer;
    private DecimalFormat df2 = new DecimalFormat("#.##", new DecimalFormatSymbols(Locale.US));
    private DaaConverter converter;
    private long maxProgress = 0, progress = 0, lastProgress = 0;

    public void write(DaaReader daaReader, String outFile, DaaConverter converter) {
        this.converter = converter;
        try {
            writer = new BufferedWriter(new FileWriter(outFile));
            writeHeader();
            ConcurrentHashMap<String, ArrayList<DaaHit>> read2hits = daaReader.getReadId2Hits();
            maxProgress = read2hits.keySet().size();
            for (String readName : read2hits.keySet()) {
                ArrayList<DaaHit> daaHits = read2hits.get(readName);
                writeReadInfo(readName, daaHits.get(0).getTotalQueryLength());
                for (DaaHit daaHit : daaHits)
                    writeSubjectAlignment(daaHit, daaReader);
                reportProgress(1, converter);
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeHeader() throws IOException {
        writer.write("BLASTP\n");
        writer.write("ELLA - Enhanced Local Aligner\n");
        writer.newLine();
    }

    public void writeReadInfo(String readName, int totalQueryLength) throws IOException {
        writer.write("Query=" + readName + "\n");
        writer.newLine();
        writer.write("Length=" + totalQueryLength + "\n");
        writer.newLine();
    }

    public void writeSubjectAlignment(DaaHit daaHit, DaaReader daaReader) throws IOException {

        DaaHeader daaHeader = daaReader.getDaaHeader();
        ScoringMatrix scoringMatrix = new ScoringMatrix(daaHeader.getScoreMatrix(), daaHeader.getGapOpen(), daaHeader.getGapExtend());
        AlignmentStatistics stats = new AlignmentStatistics(daaHit.getAlignment()[0], daaHit.getAlignment()[1], scoringMatrix);

        // writing subject header
        writer.write(">" + daaHit.getReferenceName() + "\n");
        writer.write("Length=" + daaHit.getTotalRefLength() + "\n");
        writer.newLine();

        // writing alignment scores
        writer.write("Score = " + df2.format(daaHit.getBitScore()) + " bits (" + daaHit.getRawScore() + "),\t");
        writer.write("Expected = " + df2.format(converter.getEValue(daaHit.getRawScore(), daaHit.getTotalQueryLength())) + ",\t");
        writer.newLine();
        writer.write("Identities = " + stats.getNidents() + "/" + stats.getAliLength() + " (" + df2.format(stats.getPidents()) + "%),\t");
        writer.write("Positives = " + stats.getNpositives() + "/" + stats.getAliLength() + " (" + df2.format(stats.getPpos()) + "%),\t");
        writer.write("Gaps = " + stats.getNgaps() + "/" + stats.getAliLength() + " (" + df2.format(stats.getPgaps()) + "%),\n");
        writer.write("Frame = " + daaHit.getFrame() + "\n");

        // writing alignment
        String ali1 = daaHit.getAlignment()[0], ali2 = daaHit.getAlignment()[1];
        int queryPos = daaHit.getQueryStart() + 1, refPos = daaHit.getRefStart() + 1;
        for (int l = 0; l < stats.getAliLength(); l += 60) {
            int r = Math.min(ali1.length(), l + 60);

            String sub1 = ali1.substring(l, r), sub2 = ali2.substring(l, r);
            String match = getMatchString(sub1, sub2);
            writer.newLine();

            writer.write("Query ");
            writer.write(String.format("%7d ", queryPos));
            writer.write(sub1 + " ");
            int queryLen = getSequenceLength(sub1, true);
            queryPos = daaHit.getFrame() > 0 ? queryPos + queryLen - 1 : queryPos - queryLen + 1;
            writer.write(queryPos + "");
            queryPos = daaHit.getFrame() > 0 ? queryPos + 1 : queryPos - 1;
            writer.newLine();

            writer.write("      ");
            writer.write("        ");
            writer.write(match);
            writer.newLine();

            writer.write("Sbjct ");
            writer.write(String.format("%7d ", refPos));
            writer.write(sub2 + " ");
            refPos += getSequenceLength(sub2, false) - 1;
            writer.write(refPos + "");
            refPos++;
            writer.newLine();

        }

        writer.newLine();
    }

    public void reportProgress(int delta, DaaConverter converter) {
        progress += delta;
        int p = ((int) ((((double) progress / (double) maxProgress)) * 100) / 10) * 10;
        if (p > lastProgress && p < 100) {
            lastProgress = p;
            converter.reportProgress("Converting file... ("+p+"%)");
        }
    }

    private int getSequenceLength(String aliSeq, boolean isQuery) {
        int len = 0;
        for (int i = 0; i < aliSeq.length(); i++) {
            char c = aliSeq.charAt(i);
            if (Character.isLetter(c))
                len += isQuery ? 3 : 1;
            else if (c == '\\')
                len += 1;
            else if (c == '/')
                len -= 1;
        }
        return len;
    }

    private String getMatchString(String s1, String s2) {
        StringBuilder s = new StringBuilder(s1.length());
        for (int i = 0; i < s1.length(); i++) {
            char c1 = s1.charAt(i), c2 = s2.charAt(i);
            char c = c1 != c2 ? ' ' : c1;
            s.append(c);
        }
        return s.toString();
    }

}
