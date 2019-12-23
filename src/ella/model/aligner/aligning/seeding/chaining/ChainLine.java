package ella.model.aligner.aligning.seeding.chaining;

import ella.model.aligner.aligning.seeding.Seed;

public class ChainLine {

	private StringBuffer log;
	private double m;
	private double b;

	public ChainLine(Seed[] seeds) {
		this.m = 1. / 3.;
		initLine(seeds);
	}

	private void initLine(Seed[] seeds) {

		double sumQ = 0, sumR = 0;
		for (Seed s : seeds) {
			int q = s.getQueryStart() - 1;
			int r = s.getRefStart();
			sumQ += q;
			sumR += r;
		}

		double meanQ = sumQ / (double) seeds.length;
		double meanR = sumR / (double) seeds.length;

		b = meanR - m * meanQ;

	}

	public double getDistance(Seed s) {

		int qStart = s.getQueryStart();

		double q = (double) qStart;
		double r = (double) s.getRefStart();
		double d = Math.abs(m * q - r + b) / Math.sqrt(Math.pow(m, 2));

		return d;

	}

	public double getM() {
		return m;
	}

	public double getB() {
		return b;
	}

	public StringBuffer getLog() {
		return log;
	}

}
