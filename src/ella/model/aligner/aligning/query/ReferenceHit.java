package ella.model.aligner.aligning.query;

import ella.model.aligner.aligning.seeding.chaining.SeedChain;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.MyStrand;

public class ReferenceHit {

	private MyStrand strand;
	private SeedChain seedChain;
	private long textReferenceStart;
	private int textReferenceLength;

	public ReferenceHit(MyStrand strand, long textReferenceStart, int textReferenceLength, SeedChain seedChain) {
		this.strand = strand;
		this.textReferenceStart = textReferenceStart;
		this.textReferenceLength = textReferenceLength;
		this.seedChain = seedChain;
	}

	public ReferenceHit(MyStrand strand, Object[] info, SeedChain seedChain) {
		this.strand = strand;
		this.textReferenceStart = (long) info[0];
		this.textReferenceLength = (int) info[1];
		this.seedChain = seedChain;
	}

	public ReferenceHit extractSubHit(int i) {
		SeedChain subChain = seedChain.extractSubChain(i);
		return new ReferenceHit(strand, textReferenceStart, textReferenceLength, subChain);
	}

	public double getCoverage() {
		return new Double(seedChain.getChainLength()) / new Double(textReferenceLength);
	}

	public String getAccession(IndexText text) {
		return text.getProteinAccession(textReferenceStart);
	}

	public int proteinLocationIndex(IndexText text) {
		return text.getProteinLocationIndex(textReferenceStart);
	}

	public String getSequence(IndexText text) {
		return text.getProteinSequencen(textReferenceStart);
	}

	public MyStrand getStrand() {
		return strand;
	}

	public SeedChain getSeedChain() {
		return seedChain;
	}

	public long getTextReferenceStart() {
		return textReferenceStart;
	}

	public int getTextReferenceLength() {
		return textReferenceLength;
	}

	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append("> " + textReferenceStart + " " + (textReferenceStart + textReferenceLength) + " " + strand + "\n");
		buf.append(seedChain.toString());
		return buf.toString();
	}

}
