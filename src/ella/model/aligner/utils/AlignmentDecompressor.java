package ella.model.aligner.utils;

import java.util.ArrayList;

import ella.model.io.MyParameters;

public class AlignmentDecompressor {

	public static String[] getAlignment(String queryDNA, int queryStart, ArrayList<Byte> editOperations) {

		StringBuffer positives = new StringBuffer();
		StringBuffer query = new StringBuffer();
		StringBuffer ref = new StringBuffer();
		int pos = queryStart;
		for (Byte opByte : editOperations) {
			int op = opByte & 0xFF;
			switch (op >>> 6) {
			case (0): // handling match
				for (int i = 0; i < (op & 63); i++) {
					char aa = CodonTranslator.translateCodon(queryDNA.substring(pos, pos + 3));
					query.append(aa);
					ref.append(aa);
					positives.append(aa);
					pos += 3;
				}
				break;
			case (1): // handling insertion
				for (int i = 0; i < (op & 63); i++) {
					char aa = CodonTranslator.translateCodon(queryDNA.substring(pos, pos + 3));
					query.append(aa);
					ref.append('-');
					positives.append(' ');
					pos += 3;
				}
				break;
			case (2): // handling deletion
				char c = Alphabet.getCharacter(op & 63);
				query.append('-');
				ref.append(c);
				positives.append(' ');
				break;
			case (3): // handling substitution or frameshift
				c = Alphabet.getCharacter(op & 63);
				if (c == '/' || c == '\\') {
					query.append(c);
					ref.append('-');
					positives.append(' ');
					pos += c == '/' ? -1 : 1;
				} else {
					char aa = CodonTranslator.translateCodon(queryDNA.substring(pos, pos + 3));
					ref.append(c);
					query.append(aa);
					positives.append(MyParameters.SCORING_MATRIX.getScore(c, aa) > 0 ? '|' : ' ');
					pos += 3;
				}
				break;
			}
		}

		String[] result = { query.toString(), positives.toString(), ref.toString() };
		return result;
	}

	public static boolean isOptimalAlignment(String queryDNA, int queryStart, ArrayList<Byte> editOperations) {

		int maxScore = 0, score = 0, pos = queryStart;
		boolean isGapOpen = false;
		for (int o = 0; o < editOperations.size(); o++) {
			Byte opByte = editOperations.get(o);
			int op = opByte & 0xFF;
			switch (op >>> 6) {
			case (0): // handling match
				isGapOpen = false;
				for (int i = 0; i < (op & 63); i++) {
					char aa = CodonTranslator.translateCodon(queryDNA.substring(pos, pos + 3));
					score += MyParameters.SCORING_MATRIX.getScore(aa, aa);
					pos += 3;
				}
				break;
			case (1): // handling insertion
				for (int i = 0; i < (op & 63); i++) {
					score += isGapOpen ? MyParameters.SCORING_MATRIX.getGapExtend() : MyParameters.SCORING_MATRIX.getGapOpen();
					isGapOpen = true;
					pos += 3;
				}
				break;
			case (2): // handling deletion
				score += isGapOpen ? MyParameters.SCORING_MATRIX.getGapExtend() : MyParameters.SCORING_MATRIX.getGapOpen();
				isGapOpen = true;
				break;
			case (3): // handling substitution or frameshift
				isGapOpen = false;
				char c = Alphabet.getCharacter(op & 63);
				if (c == '/' || c == '\\') {
					score += MyParameters.FRAMESHIFT_PENALTY;
					pos += c == '/' ? -1 : 1;
				} else {
					char aa = CodonTranslator.translateCodon(queryDNA.substring(pos, pos + 3));
					score += MyParameters.SCORING_MATRIX.getScore(aa, c);
					pos += 3;
				}
				break;
			}

			maxScore = maxScore < score ? score : maxScore;
			if (score < maxScore - MyParameters.X_BESTSCORE_DROP) {
				return false;
			}

		}
		return true;
	}

}
