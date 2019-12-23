package ella.model.io.converting;

import ella.model.aligner.utils.ScoringMatrix;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class BlastTabFile {

    private DecimalFormat df2 = new DecimalFormat("#.##");
    private long maxProgress = 0, progress = 0, lastProgress = 0;

    public void write(DaaReader daaReader, String outFile, int cores, boolean writeHeaderLine, DaaConverter converter) {

        try {

            BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
            String headerLine = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore \n";
            if (writeHeaderLine)
                writer.write(headerLine);
            DaaHeader daaHeader = daaReader.getDaaHeader();
            ScoringMatrix scoringMatrix = new ScoringMatrix(daaHeader.getScoreMatrix(), daaHeader.getGapOpen(), daaHeader.getGapExtend());

            ArrayList<DaaHit> daaHits = daaReader.parseAllHits(cores);
            maxProgress = daaHits.size();
            for (DaaHit hit : daaHits) {

                String ali1 = hit.getAlignment()[0];
                String ali2 = hit.getAlignment()[1];

                String qseqid = hit.getQueryName(); // qseqid
                int qlen = hit.getTotalQueryLength(); // qlen

                String sseqid = hit.getReferenceName(); // sseqid
                int slen = hit.getTotalRefLength(); //slen

                int qstart = hit.getQueryStart() + 1; // qstart
                int qend = hit.getFrame() < 0 ? hit.getQueryStart() - hit.getQueryLength() + 2 : hit.getQueryStart() + hit.getQueryLength();

                int sstart = hit.getRefStart() + 1; // sstart
                int send = hit.getRefStart() + hit.getRefLength(); // send

                String evalue = df2.format(converter.getEValue(hit.getRawScore(), hit.getTotalQueryLength())); // evalue
                String bitscore = df2.format(converter.getBitScore(hit.getRawScore())); // bitscore
                int score = hit.getRawScore(); // score
                int qframe = hit.getFrame(); // qframe

                AlignmentStatistics stats = new AlignmentStatistics(ali1, ali2, scoringMatrix);
                String pident = df2.format(stats.getPidents());

                String lineToWrite = qseqid + "\t" + sseqid + "\t" + pident + "\t" + stats.getAliLength() + "\t" + stats.getMismatches() + "\t" + stats.getNgops() + "\t" + qstart + "\t" + qend + "\t" + sstart + "\t" + send + "\t" + evalue + "\t" + bitscore + "\n";

                writer.write(lineToWrite);
                reportProgress(1, converter);

            }

            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void reportProgress(int delta, DaaConverter converter) {
        progress += delta;
        int p = ((int) ((((double) progress / (double) maxProgress)) * 100) / 10) * 10;
        if (p > lastProgress && p < 100) {
            lastProgress = p;
            converter.reportProgress("Converting file... (" + p + "%)");
        }
    }

}
