package ella.model.aligner.aligning.seeding.chaining;

import java.util.ArrayList;

import ella.model.aligner.aligning.seeding.Seed;
import ella.model.aligner.utils.MyStrand;

public class SeedChain {

	private ArrayList<Seed> seeds;
	private MyStrand strand;

	public SeedChain(ArrayList<Seed> seeds) {
		this.seeds = seeds;
		if (!seeds.isEmpty())
			strand = seeds.get(0).getStrand();
	}

	public SeedChain extractSubChain(int pos) {
		ArrayList<Seed> subSeeds = new ArrayList<Seed>();
		for (int i = pos; i < seeds.size(); i++)
			subSeeds.add(seeds.get(i));
		return new SeedChain(subSeeds);
	}

	public int getChainLength() {
		int length = 0;
		for (Seed s : seeds)
			length += s.getRefLength();
		return length;
	}

	public int getSize() {
		return seeds.size();
	}

	public ArrayList<Seed> getSeeds() {
		return seeds;
	}

	public MyStrand getStrand() {
		return strand;
	}

	public boolean isEmpty() {
		return seeds.isEmpty();
	}

	public String toString() {
		StringBuffer buf = new StringBuffer();
		for (Seed seed : seeds) {
			int refStart = seed.getRefStart();
			int refEnd = seed.getRefStart() + seed.getRefLength();
			int queryStart = seed.getQueryStart();
			int queryEnd = seed.getQueryStart() + seed.getQueryLength();
			buf.append("[" + refStart + "," + refEnd + "]\t" + "[" + queryStart + "," + queryEnd + "]\t" + seed.getFrame() + "\n");
		}
		return buf.toString();
	}

}
