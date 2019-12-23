package ella.model.io.converting;

import ella.model.aligner.utils.ScoringMatrix;

class AlignmentStatistics {

    private String ali1, ali2;
    private ScoringMatrix scoringMatrix;
    private int npositives = 0, nidents = 0, mismatches = 0, aliLength = 0;
    private int ngops = 0, ngaps = 0;
    private double pidents, ppos, pgaps;

    public AlignmentStatistics(String ali1, String ali2, ScoringMatrix scoringMatrix) {
        this.ali1 = ali1;
        this.ali2 = ali2;
        this.scoringMatrix = scoringMatrix;
        compute();
    }

    private void compute() {
        aliLength = ali1.length(); // length
        npositives = nidents = mismatches = 0; // mismatch
        ngops = ngaps = 0;
        boolean insideGap = false;
        for (int i = 0; i < aliLength; i++) {
            char c1 = ali1.charAt(i), c2 = ali2.charAt(i);
            if (Character.isLetter(c1) && Character.isLetter(c2)) {
                if (scoringMatrix.getScore(c1, c2) > 0)
                    npositives++;
                if (ali1.charAt(i) == ali2.charAt(i)) {
                    nidents++;
                } else {
                    mismatches++;
                }
                insideGap = false;
            } else if (Character.isLetter(c1) || Character.isLetter(c2)) {
                if (!insideGap) {
                    ngaps++;
                    ngops++;
                    insideGap = true;
                } else
                    ngaps++;
            } else
                insideGap = false;
        }

        pidents = (((double) nidents * 100.) / (double) aliLength);// pident
        ppos = (((double) npositives * 100.) / (double) aliLength);// pident
        pgaps = (((double) ngaps * 100.) / (double) aliLength);

    }

    public double getPgaps() {
        return pgaps;
    }

    public int getAliLength() {
        return aliLength;
    }

    public int getNpositives() {
        return npositives;
    }

    public int getNidents() {
        return nidents;
    }

    public int getMismatches() {
        return mismatches;
    }

    public int getNgops() {
        return ngops;
    }

    public int getNgaps() {
        return ngaps;
    }

    public double getPidents() {
        return pidents;
    }

    public double getPpos() {
        return ppos;
    }
}
