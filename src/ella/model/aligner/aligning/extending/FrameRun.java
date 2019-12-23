package ella.model.aligner.aligning.extending;

import java.util.ArrayList;

import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.aligning.query.ReferenceHit;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.MyStrand;

public class FrameRun {

	private ArrayList<FrameHit> frameHits;
	private FrameShiftHit frameShiftHit;
	private MyStrand strand;
	private int refRunLength, proteinLocationIndex;
	private long textReferenceStart;
	private int textReferenceLength;
	private QueryContainer qC;

	public FrameRun(ArrayList<FrameHit> frameHits, ReferenceHit refHit, QueryContainer qC, IndexText text) {
		this.frameHits = frameHits;
		this.strand = refHit.getStrand();
		this.proteinLocationIndex = refHit.proteinLocationIndex(text);
		this.qC = qC;
		this.textReferenceStart = refHit.getTextReferenceStart();
		this.textReferenceLength = refHit.getTextReferenceLength();
		FrameHit firstHit = frameHits.get(0);
		FrameHit lastHit = frameHits.get(frameHits.size() - 1);
		this.refRunLength = lastHit.getRefEnd() - firstHit.getRefStart();
		this.frameShiftHit = new FrameShiftHit(frameHits, textReferenceStart);
	}

	public ArrayList<FrameHit> getFrameHits() {
		return frameHits;
	}

	public MyStrand getStrand() {
		return strand;
	}

	public int getRefCoverage() {
		return (int) Math.floor(((double) refRunLength / (double) textReferenceLength) * 100.);
	}

	public String getRefAccession(IndexText text) {
		return text.getProteinAccession(textReferenceStart);
	}

	public String getRefSequence(IndexText text) {
		return text.getProteinSequencen(textReferenceStart);
	}

	public int getRefRunLength() {
		return refRunLength;
	}

	public String getReadName() {
		return qC.getId().toString();
	}

	public byte[] getPackedSequence() {
		return qC.getPackedSequence();
	}

	public int getTotalQueryLenth() {
		return qC.getTotalQueryLength();
	}

	public FrameShiftHit getFrameShiftHit() {
		return frameShiftHit;
	}

	public int getProteinLocationIndex() {
		return proteinLocationIndex;
	}

	public double getBitScore() {
		return frameShiftHit.getBitScore();
	}

	public double getEValue() {
		return frameShiftHit.cmpEValue(qC.getTotalQueryLength());
	}

	public double getEG2Value() {
		return frameShiftHit.cmpEG2Value();
	}

	public double getRawScore() {
		return frameShiftHit.getRawScore();
	}

	public boolean checkConsistency(String info) {
		for (int i = 1; i < frameHits.size(); i++) {
			FrameHit f1 = frameHits.get(i - 1);
			FrameHit f2 = frameHits.get(i);
			if (f1.getRefEnd() != f2.getRefStart() || Math.abs(f1.getQueryEnd() - f2.getQueryStart()) > 1) {
				System.out.println("!!! FRAME RUN INCONSISTENT !!! ");
				System.out.println(">" + info);
				System.out.println(f1.toString());
				System.out.println(f2.toString());
				return false;
			}
		}
		return true;
	}

	public String toString(IndexText text) {
		StringBuffer buf = new StringBuffer();
		buf.append(">" + getRefAccession(text) + " Cov: " + getRefCoverage() + " RawScore: " + getRawScore() + " BitScore: " + getBitScore()
				+ " EValue: " + getEValue() + " EG2Value: " + getEG2Value() + "\n");
		for (String s : frameShiftHit.getAlignment(qC.getStrandSpecificDNA(strand)))
			buf.append(s + "\n");
		return buf.toString();
	}

}
