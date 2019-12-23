package ella.model.aligner.aligning.seeding;

import java.util.ArrayList;

public class DiagonalTable {

    private final int n = 64;
    private ArrayList<long[]>[] coordinateBins = new ArrayList[n];

    public boolean coversCoordinate(int x, long y) {
        long diagonal = x - y;
        int b = getBucketID(diagonal);
        if (coordinateBins[b] != null) {
            for (long[] entry : coordinateBins[b]) {
                if (entry[0] > x && entry[1] == diagonal) {
                    return true;
                }
            }
        }
        return false;
    }

    public void addCoordinate(int x, long y) {
        long diagonal = x - y;
        int b = getBucketID(diagonal);
        long[] coord = {x, diagonal};
        if (coordinateBins[b] == null)
            coordinateBins[b] = new ArrayList<>();
        coordinateBins[b].add(coord);
    }

    private int getBucketID(long d) {
        int b = (int) (d - (d / n * n));
        b = b < 0 ? -b : b;
        return b;
    }

}
