package ella.model.aligner.utils.export;

import ella.model.aligner.utils.daaMerger.DaaHit;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;

public class BlastTab {

    public static void run(File out, ArrayList<DaaHit> allHits) {

        try {

            BufferedWriter writer = new BufferedWriter(new FileWriter(out));

            String header = "qseqid" + "\t" + "qlen" + "\t" + "sseqid" + "\t" + "slen" + "\t"
                    + "qstart" + "\t" + "qend" + "\t" + "sstart" + "\t" + "send" + "\t" + "evalue" + "\t" + "bitscore" + "\t" + "score" + "\t" + "length" + "\t"
                    + "pident" + "\t" + "nident" + "\t" + "mismatch" + "\t" + "positive" + "\t" + "gapopen" + "\t"
                    + "gaps" + "\t" + "ppos" + "\t" + "qframe" + "\t" + "qseq" + "\t" + "sseq" + "\n";
            writer.write(header);

            for (DaaHit hit : allHits) {

                String seq1 = hit.getAlignment()[0];
                String seq2 = hit.getAlignment()[1];

                String qseqid = hit.getQueryName(); // qseqid
                int qlen = hit.getTotalQueryLength(); // qlen

                String sseqid = hit.getReferenceName(); // sseqid
//            String sallseqid = getAccession(hit); // sallseqid
                int slen = hit.getTotalRefLength(); //slen

                int qstart = hit.getQueryStart() + 1; // qstart
                int qend;
                if (hit.getFrame() < 0) {
                    qend = hit.getQueryStart() - hit.getQueryLength() + 2; // qend
                } else {
                    qend = hit.getQueryStart() + hit.getQueryLength(); // qend
                }

                int sstart = hit.getRefStart() + 1; // sstart
                int send = hit.getRefStart() + hit.getRefLength(); // send

                String qseq = hit.getQueryDNA(); // qseq
                String sseq = seq2.replaceAll("-", ""); // sseq

                double evalue = AlignmentStatisticsHelper.getEValue(hit.getRawScore(), hit.getQueryLength()); // evalue
                double bitscore = AlignmentStatisticsHelper.getBitScore(hit.getRawScore()); // bitscore
                int score = hit.getRawScore(); // score
                int length = hit.getRefLength(); // length

                int nident = 0; // nident
                int mismatch = 0; // mismatch
                for (int i = 0; i < length; i++) {
                    if (seq1.charAt(i) != '-' && seq1.charAt(i) != '/' && seq1.charAt(i) != '\\' && seq2.charAt(i) != '-' && seq2.charAt(i) != '/' && seq2.charAt(i) != '\\') {
                        if (seq1.charAt(i) == seq2.charAt(i)) {
                            nident++;
                        } else {
                            mismatch++;
                        }
                    }
                }

                double pident = ((double) nident / (double) length) * 100; // pident

//                int positives = hit.getPositives(); // positive (for now this is not used)
                int positives = 0; // TODO

                int gapopen = 0; // gapopen
                int gaps = 0; // gaps

                boolean insideGap = false;
                for (int i = 0; i < length; i++) {
                    if (seq1.charAt(i) == '-' && !insideGap) {
                        insideGap = true;
                        gapopen++;
                    } else if (seq1.charAt(i) == '-' && insideGap) {
                        gaps++;
                    } else if (seq1.charAt(i) != '-') {
                        insideGap = false;
                    }
                }

                insideGap = false;
                for (int i = 0; i < length; i++) {
                    if (seq2.charAt(i) == '-' && !insideGap) {
                        insideGap = true;
                        gapopen++;
                        gaps++;
                    } else if (seq2.charAt(i) == '-' && insideGap) {
                        gaps++;
                    } else if (seq2.charAt(i) != '-') {
                        insideGap = false;
                    }
                }

                double ppos = ((double) positives / (double) length) * 100; // ppos

                int qframe = hit.getFrame(); // qframe

//            String btop = null; // btop

//            String staxids = null; // staxids

//            String stitle = hit.getReferenceName(); // stitle
//            String salltitles = hit.getReferenceName(); // salltitles

//            double qcovhsp = 0.0; // qcovhsp

//            String qtitle = hit.getQueryName(); // qtitle

                String lineToWrite = qseqid + "\t" + qlen + "\t" + sseqid + "\t" + slen + "\t" + qstart
                        + "\t" + qend + "\t" + sstart + "\t" + send + "\t"
                        + evalue + "\t" + bitscore + "\t" + score + "\t" + length + "\t" + pident + "\t" + nident + "\t"
                        + mismatch + "\t" + positives + "\t" + gapopen + "\t" + gaps + "\t" + ppos + "\t" + qframe + "\t" + qseq + "\t" + sseq +
                        "\n";

                writer.write(lineToWrite);

            }

            writer.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
