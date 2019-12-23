package ella.model.aligner.aligning.seeding;

import ella.model.aligner.aligning.seeding.chaining.ChainLine;
import ella.model.aligner.utils.MyFrame;
import ella.model.aligner.utils.MyStrand;

public class Seed {

	private MyFrame frame;
	private long indexStart;
	private int queryStart, refStart, refLength;
	private int score;

	public Seed(MyFrame frame, int queryStart, int refStart, long indexStart, int refLength, int score) {
		this.frame = frame;
		this.queryStart = queryStart;
		this.indexStart = indexStart;
		this.refStart = refStart;
		this.refLength = refLength;
		this.score = score;
	}

	public void increaseQueryStart(int offset) {
		offset = Math.min(offset, getQueryLength());
		queryStart += offset;
		refStart += (offset / 3);
		refLength -= (offset / 3);
	}

	public boolean contains(Seed s) {
		return (s.getQueryStart() >= queryStart && s.getQueryEnd() <= getQueryEnd() && s.getIndexStart() >= indexStart
				&& s.getIndexEnd() <= getIndexEnd());
	}

	public int getYCoordinate(int referenceLength) {
		int y = queryStart - 3 * refStart;
		y += referenceLength * 3;
		return y;
	}

	public int getBinIndex(double epsilon, int referenceLength) {
		double y = getYCoordinate(referenceLength);
		int index = (int) Math.floor(y / epsilon);
		return index;
	}

	public double getDistToLine(ChainLine hitLine) {
		return hitLine.getDistance(this);
	}

	public MyStrand getStrand() {
		return frame.getStrand();
	}

	public MyFrame getFrame() {
		return frame;
	}

	public int getRefStart() {
		return refStart;
	}

	public long getIndexStart() {
		return indexStart;
	}

	public long getIndexEnd() {
		return indexStart + refLength;
	}

	public int getRefEnd() {
		return refStart + refLength;
	}

	public int getRefLength() {
		return refLength;
	}

	public int getQueryStart() {
		return queryStart;
	}

	public int getQueryEnd() {
		return queryStart + refLength * 3;
	}

	public int getQueryLength() {
		return refLength * 3;
	}

	public int getScore() {
		return score;
	}

	public String toString() {
		return "[" + indexStart + "," + getIndexEnd() + "]\t" + "[" + refStart + "," + getRefEnd() + "]\t" + "[" + queryStart + "," + getQueryEnd()
				+ "]\t" + frame + " " + score;
	}

}
