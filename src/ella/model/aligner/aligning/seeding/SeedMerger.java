package ella.model.aligner.aligning.seeding;

import java.util.ArrayList;
import java.util.Comparator;

public class SeedMerger {

	public static ArrayList<Seed> run(Seed[] seeds) {

		if (seeds.length == 0)
			return new ArrayList<Seed>();

		// merging seeds of the same frame
		ArrayList<Seed> mergedSeeds = new ArrayList<Seed>();
		ArrayList<Seed> seedSet = new ArrayList<Seed>();
		seedSet.add(seeds[0]);
		for (int i = 1; i < seeds.length; i++) {
			Seed seed2 = seeds[i];
			if (overlaps(seedSet.get(seedSet.size() - 1), seed2)) {
				seedSet.add(seed2);
			} else {
				mergedSeeds.addAll(mergeSeeds(seedSet));
				seedSet.clear();
				seedSet.add(seed2);
			}
		}
		mergedSeeds.addAll(mergeSeeds(seedSet));

		// removing encapsulated seeds
		ArrayList<Seed> toDelete = new ArrayList<Seed>();
		do {
			toDelete = new ArrayList<Seed>();
			for (int i = 1; i < mergedSeeds.size(); i++) {
				Seed seed1 = mergedSeeds.get(i - 1);
				Seed seed2 = mergedSeeds.get(i);
				if (seed1.contains(seed2))
					toDelete.add(seed2);
			}
			for (Seed s : toDelete)
				mergedSeeds.remove(s);
		} while (!toDelete.isEmpty());

		// resolving seed overlaps amongst different frames
		for (int i = 1; i < mergedSeeds.size(); i++) {
			Seed seed1 = mergedSeeds.get(i - 1);
			Seed seed2 = mergedSeeds.get(i);
			if (seed2.getQueryStart() < seed1.getQueryEnd()) {
				int offset = seed1.getQueryEnd() - seed2.getQueryStart();
				while (offset % 3 != 0)
					offset++;
				seed2.increaseQueryStart(offset);
				if (seed2.getQueryStart() == seed2.getQueryEnd())
					toDelete.add(seed2);
			}
		}
		for (Seed s : toDelete)
			mergedSeeds.remove(s);

		return mergedSeeds;
	}

	// merging overlapping seeds of the same frame should actually never happen
	private static ArrayList<Seed> mergeSeeds(ArrayList<Seed> seedSet) {
		if (seedSet.size() < 2)
			return seedSet;

		Seed seed1 = seedSet.get(0);
		Seed seed2 = seedSet.get(seedSet.size() - 1);
		int length = seed2.getRefEnd() - seed1.getRefStart();
		// TODO: there is a better way to do that instead of simply choosing the max
		int score = seedSet.stream().max(Comparator.comparing(Seed::getScore)).get().getScore();

		Seed seed = new Seed(seed1.getFrame(), seed1.getQueryStart(), seed1.getRefStart(), seed1.getIndexStart(), length, score);
		seedSet.clear();
		seedSet.add(seed);

		return seedSet;
	}

	private static boolean overlaps(Seed seed1, Seed seed2) {
		if (seed1.getFrame() != seed2.getFrame())
			return false;
		int qDiff = seed2.getQueryStart() - seed1.getQueryStart();
		int rDiff = seed2.getRefStart() - seed1.getRefStart();
		return (qDiff == rDiff * 3 && qDiff < seed1.getQueryLength());
	}

	public static class SeedSorter implements Comparator<Seed> {

		@Override
		public int compare(Seed s1, Seed s2) {
			return Integer.compare(s1.getQueryStart(), s2.getQueryStart());
		}

	}

}
