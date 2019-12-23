package ella.model.aligner.utils;

import java.io.*;
import java.util.ArrayList;

public class ScoringMatrix {

    private String type;
    private int gapOpen, gapExtend;
    private ArrayList<String> alphabet;
    private int[][] matrix;

    public ScoringMatrix(String type, int gapOpen, int gapExtend) {
        this.type = type;
        this.gapOpen = gapOpen;
        this.gapExtend = gapExtend;
        parseMatrix();
    }

    public int[] cmpAlignmentScores(String s1, String s2) {
        boolean isGapOpen = false;
        int[] scores = new int[s1.length()];
        for (int i = 0; i < s1.length(); i++) {

            char a = s1.charAt(i);
            char b = s2.charAt(i);
            int score = -gapOpen - gapExtend;
            if (a == '-' || b == '-') {
                if (isGapOpen)
                    score = -gapExtend;
                isGapOpen = true;
            } else {
                score = getScore(a, b);
                isGapOpen = false;
            }

            scores[i] = score;

        }
        return scores;
    }

    public int cmpAlignmentScore(String s1, String s2) {
        boolean isGapOpen = false;
        int score = 0;
        for (int i = 0; i < s1.length(); i++) {
            char a = s1.charAt(i);
            char b = s2.charAt(i);
            int s = -gapOpen - gapExtend;
            if (a == '-' || b == '-') {
                if (isGapOpen)
                    s = -gapExtend;
                isGapOpen = true;
            } else {
                s = getScore(a, b);
                isGapOpen = false;
            }

            score += s;
        }
        return score;
    }

    public int getScore(char a, char b) {
        int i = getIndex(a);
        int j = getIndex(b);
        return matrix[i][j];
    }

    private int getIndex(char c) {
        int i = (int) c;
        if (i >= 97 && i <= 122)
            return i - 65 - 32;
        if (i >= 65 && i <= 90)
            return i - 65;
        return 26;
    }

    private int parseMatrix() {

        try {

            InputStream is = this.getClass().getResourceAsStream("/NCBI_ScoringMatrices/" + type);
            BufferedReader buf = new BufferedReader(new InputStreamReader(is));

            alphabet = new ArrayList<String>();
            for (int i = 0; i < 2; i++) {
                String l = buf.readLine();
                ArrayList<String> entries = splitLine(l);
                alphabet.addAll(entries);
            }

            matrix = new int[27][27];
            int i = 0;
            String l;
            ArrayList<String> buffer = new ArrayList<String>();
            while ((l = buf.readLine()) != null) {
                if (!l.isEmpty()) {
                    buffer.addAll(splitLine(l));
                    if (buffer.size() == alphabet.size()) {
                        int j = 0;
                        int row = getIndex(alphabet.get(i).charAt(0));
                        for (String e : buffer) {
                            int val = Integer.parseInt(e);
                            int col = getIndex(alphabet.get(j).charAt(0));
                            matrix[row][col] = val;
                            j++;
                        }
                        buffer.clear();
                        i++;
                    }
                }

            }
            buf.close();

            return 1;

        } catch (Exception e) {
            e.printStackTrace();
        }

        return 0;
    }

    private ArrayList<String> splitLine(String l) {
        String[] entries = l.split("\\s+");
        ArrayList<String> splitResult = new ArrayList<String>();
        for (String e1 : entries) {
            for (String e2 : e1.split("\\,")) {
                if (!e2.isEmpty() && !e2.startsWith("/*") && !e2.startsWith("*/")) {
                    splitResult.add(e2.replaceAll("/,", ""));
                }
            }
        }
        return splitResult;
    }

    public int getGapOpen() {
        return gapOpen;
    }

    public int getGapExtend() {
        return gapExtend;
    }

    public String getType() {
        return type;
    }

    public int[][] getMatrix() {
        return matrix;
    }

    public void print() {

        String out = type + " Matrix:\n";
        out = "\t";
        for (String c : alphabet)
            out = out.concat(c + "\t");
        out = out.concat("\n");

        for (int i = 0; i < alphabet.size(); i++) {
            out = out.concat(alphabet.get(i) + "\t");
            for (int j = 0; j < alphabet.size(); j++) {
                int row = getIndex(alphabet.get(i).charAt(0));
                int col = getIndex(alphabet.get(j).charAt(0));
                out = out.concat(matrix[row][col] + "\t");
            }
            out = out.concat("\n");
        }
        System.out.println(out);

    }

}
