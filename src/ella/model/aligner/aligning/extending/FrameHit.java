package ella.model.aligner.aligning.extending;

import java.util.ArrayList;

public class FrameHit {

	private int frame;
	private int score;
	private int queryStart, refStart, queryEnd, refEnd;
	private ArrayList<Byte> editOperations;

	public FrameHit(int frame, int score, int queryStart, int refStart, int queryLength, int refLength, ArrayList<Byte> editOperations) {
		super();
		this.frame = frame;
		this.score = score;
		this.queryStart = queryStart;
		this.refStart = refStart;
		this.queryEnd = queryStart + queryLength;
		this.refEnd = refStart + refLength;
		this.editOperations = editOperations;
	}

	public int getFrame() {
		return frame;
	}

	public int getScore() {
		return score;
	}

	public int getQueryStart() {
		return queryStart;
	}

	public int getRefStart() {
		return refStart;
	}

	public int getQueryEnd() {
		return queryEnd;
	}

	public int getRefEnd() {
		return refEnd;
	}

	public ArrayList<Byte> getEditOperations() {
		return editOperations;
	}

	public String toString() {
		return "[" + refStart + "," + refEnd + "]\t" + "[" + queryStart + "," + queryEnd + "]\t" + frame + "\t" + score;
	}

	public FrameHit addFrameHit(FrameHit f2) {
		ArrayList<Byte> mergedEditOperations = new ArrayList<Byte>();
		mergedEditOperations.addAll(editOperations);
		mergedEditOperations.addAll(f2.getEditOperations());
		return new FrameHit(frame, score + f2.getScore(), queryStart, refStart, f2.getQueryEnd() - queryStart, f2.getRefEnd() - refStart,
				mergedEditOperations);
	}

}
