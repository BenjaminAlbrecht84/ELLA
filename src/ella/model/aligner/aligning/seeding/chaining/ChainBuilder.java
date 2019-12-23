package ella.model.aligner.aligning.seeding.chaining;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.aligning.query.ReferenceHit;
import ella.model.aligner.aligning.seeding.Seed;
import ella.model.aligner.aligning.seeding.SeedMerger;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.MyStrand;
import ella.model.aligner.utils.streams.MySeedStream;
import ella.model.aligner.utils.wrapper.ObjectArrayWrapper;
import ella.model.io.MyParameters;

public class ChainBuilder {

    private final SeedRadixSort seedRadixSort = new SeedRadixSort();

    public void run(QueryContainer qC, IndexText text) {

        // assessing best seed-set in plus strand
        HashMap<ObjectArrayWrapper, SeedChain> seedMap = new HashMap<ObjectArrayWrapper, SeedChain>();
        for (ObjectArrayWrapper key : qC.getPlusSeedMap().keySet()) {
            SeedChain plusChain = chainSeeds(qC.getPlusSeedMap().get(key).toArray(), MyStrand.FWD, (int) key.getData()[1],
                    qC.getTotalQueryLength());
            if (!plusChain.isEmpty())
                seedMap.put(key, plusChain);
        }

		// assessing best seed-set in minus strands
        for (ObjectArrayWrapper key : qC.getMinusSeedMap().keySet()) {
            SeedChain minusChain = chainSeeds(qC.getMinusSeedMap().get(key).toArray(), MyStrand.REV, (int) key.getData()[1],
                    qC.getTotalQueryLength());
            double minusLength = minusChain.getChainLength();
            double plusLength = seedMap.containsKey(key) ? seedMap.get(key).getChainLength() : 0;
            if (minusLength > plusLength)
                seedMap.put(key, minusChain);
        }

        // storing best seed-set result
        for (ObjectArrayWrapper key : seedMap.keySet()) {
            SeedChain seedChain = seedMap.get(key);
            if (!seedMap.get(key).isEmpty()) { // && key.getData()[1] > 50) {
                ReferenceHit rH = new ReferenceHit(seedChain.getStrand(), key.getData(), seedChain);
                qC.addRefHit(rH);
            }

        }

    }

    private void plotSeeds(ArrayList<Seed> seeds) {
        StringBuilder x = new StringBuilder("x <- c(0");
        StringBuilder y = new StringBuilder("y <- c(0");
        for (Seed s : seeds) {
            x.append("," + s.getQueryStart());
            y.append("," + s.getRefStart());
        }
        System.out.println("\n" + x + ")");
        System.out.println(y + ")");
        System.out.println("plot(x,y)");
    }

    public SeedChain chainSeeds(Seed[] inputSeeds, MyStrand strand, int referenceLength, int queryLength) {

        ArrayList<Seed[]> allSeeds = new ArrayList<>();
        if (MyParameters.FRAMESHIFT_PENALTY != null)
            allSeeds.add(inputSeeds);
        else {
            HashMap<Integer, ArrayList<Seed>> frame2seeds = new HashMap<>();
            for (Seed seed : inputSeeds) {
                int frame = seed.getFrame().getNumericalID();
                frame2seeds.putIfAbsent(frame, new ArrayList<>());
                frame2seeds.get(frame).add(seed);
            }
            for (int frame : frame2seeds.keySet()) {
                Seed[] seedArray = frame2seeds.get(frame).toArray(new Seed[frame2seeds.get(frame).size()]);
                allSeeds.add(seedArray);
            }
        }

        SeedChain bestSeedChain = null;
        for (Seed[] seeds : allSeeds) {

            double epsilon = Math.max(0.06 * (double) referenceLength, 10.0);

            seeds = seedRadixSort.sort(seeds, referenceLength);

            // computing set of binned seeds
            Seed[] binnedSeeds = seeds; // binSeeds(seeds, referenceLength, epsilon);
            if (binnedSeeds.length == 0)
                return new SeedChain(SeedMerger.run(binnedSeeds));

            // filtering binned seeds
            int[] s = {0, 0, 0};
            int[] sMax = {0, 1, 0};
            while (s[1] < binnedSeeds.length) {

                // extending seed section
                while (s[1] < binnedSeeds.length && dist(binnedSeeds[s[0]], binnedSeeds[s[1]], referenceLength) <= epsilon) {
                    s[2] += binnedSeeds[s[1]].getScore();
                    s[1]++;
                }

                // storing best seed section
                if (s[1] - s[0] > sMax[1] - sMax[0] || (s[1] - s[0] == sMax[1] - sMax[0] && s[2] > sMax[2])) {
                    sMax[0] = s[0];
                    sMax[1] = s[1];
                    sMax[2] = s[2];
                }

                while (s[1] < binnedSeeds.length && dist(binnedSeeds[s[0]], binnedSeeds[s[1]], referenceLength) > epsilon) {
                    s[2] -= binnedSeeds[s[0]].getScore();
                    s[0]++;
                }

            }

            // linearizing best seed section
            Seed[] linearSeeds = new Seed[sMax[1] - sMax[0]];
            int pos = 0;
            for (int i = sMax[0]; i < sMax[1]; i++)
                linearSeeds[pos++] = binnedSeeds[i];
            Arrays.sort(linearSeeds, new SeedXComparator());

            linearSeeds = linearizeSeedSet(linearSeeds, new ChainLine(linearSeeds));
            SeedChain seedChain = new SeedChain(SeedMerger.run(linearSeeds));
            bestSeedChain = bestSeedChain == null || seedChain.getChainLength() > bestSeedChain.getChainLength() ? seedChain : bestSeedChain;

        }

        return bestSeedChain;

    }

    private Seed[] linearizeSeedSet(Seed[] linearSeeds, ChainLine hitLine) {

        ArrayList<Seed> toRemove = new ArrayList<Seed>();
        ArrayList<Seed> toRefine = new ArrayList<Seed>();
        for (int i = 1; i < linearSeeds.length; i++) {
            Seed s1 = linearSeeds[i - 1];
            Seed s2 = linearSeeds[i];
            if (s1.getRefEnd() > s2.getRefStart() || s1.getQueryEnd() > s2.getQueryStart()) {
                if (!toRefine.contains(s1))
                    toRefine.add(s1);
                toRefine.add(s2);
            } else if (!toRefine.isEmpty()) {
                refineBlockRec(toRefine, toRemove, hitLine);
                toRefine.clear();
            }
        }
        if (!toRefine.isEmpty())
            refineBlockRec(toRefine, toRemove, hitLine);
        if (!toRemove.isEmpty()) {
            MySeedStream seedStream = new MySeedStream(linearSeeds.length);
            for (Seed s : linearSeeds) {
                if (!toRemove.contains(s))
                    seedStream.add(s);
            }
            return linearizeSeedSet(seedStream.toArray(), hitLine);
        }

        return linearSeeds;
    }

    private void refineBlockRec(ArrayList<Seed> toRefine, ArrayList<Seed> toRemove, ChainLine hitLine) {
        Seed delSeed = toRefine.stream().min(Comparator.comparing(s -> s.getScore())).get();
        toRemove.add(delSeed);
        toRefine.remove(delSeed);
        boolean isLinear = true;
        for (int i = 1; i < toRefine.size(); i++) {
            Seed s1 = toRefine.get(i - 1);
            Seed s2 = toRefine.get(i);
            if (s1.getRefStart() > s2.getRefStart()) {
                isLinear = false;
                break;
            }
        }
        if (!isLinear)
            refineBlockRec(toRefine, toRemove, hitLine);
    }

    private double dist(Seed s1, Seed s2, int referenceLength) {
        double y1 = s1.getYCoordinate(referenceLength);
        double y2 = s2.getYCoordinate(referenceLength);
        return Math.abs(y1 - y2);
    }

    public class SeedYComparator implements Comparator<Seed> {

        private int referenceLength;

        public SeedYComparator(int referenceLength) {
            this.referenceLength = referenceLength;
        }

        @Override
        public int compare(Seed s1, Seed s2) {
            double y1 = s1.getYCoordinate(referenceLength);
            double y2 = s2.getYCoordinate(referenceLength);
            return Double.compare(y1, y2);
        }

    }

    public class SeedXComparator implements Comparator<Seed> {
        @Override
        public int compare(Seed s1, Seed s2) {
            int x1 = s1.getQueryStart();
            int x2 = s2.getQueryStart();
            return Integer.compare(x1, x2);
        }
    }

}
