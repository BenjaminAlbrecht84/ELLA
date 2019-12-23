package ella.model.aligner.utils;

import java.util.LinkedList;

public class Minimizer {

	private int maxSize;

	private LinkedList<Object[]> pMerList;
	private int minScore, minIndex;

	public Minimizer(int p, int q) {
		this.maxSize = q - p + 1;
		pMerList = new LinkedList<Object[]>();
		minScore = Integer.MAX_VALUE;
	}

	public void addPMer(Object o, String p) {

		// adding new pMer and updating minimum
		Object[] pMer = { o, cmpScore(p) };
		pMerList.addLast(pMer);
		if ((int) pMer[1] <= minScore) {
			minIndex = pMerList.size() - 1;
			minScore = (int) pMer[1];
		}

		// kicking out first one and assessing new minimum if necessary
		if (pMerList.size() > maxSize) {
			if (minIndex == 0)
				findNewMinimium();
			pMerList.removeFirst();
			minIndex--;
		}

	}

	private int cmpScore(String p) {
		int score = 0;
		for (int i = 0; i < p.length(); i++)
			score += Alphabet.reduceCharacter(p.charAt(i));
		return score;
	}

	private void findNewMinimium() {
		minScore = Integer.MAX_VALUE;
		// first one will be kicked out, so it's ignored here
		for (int i = 1; i < pMerList.size(); i++) {
			Object[] pMer = pMerList.get(i);
			if ((int) pMer[1] < minScore) {
				minScore = (int) pMer[1];
				minIndex = i;
			}
		}
	}

	public Object getMinimium() {
		if (pMerList.size() <= minIndex)
			return null;
		return pMerList.get(minIndex)[0];
	}

}
