/*
 * Copyright 2017 Benjamin Albrecht
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package ella.model.aligner.utils.daaMerger;

import java.util.ArrayList;

import ella.model.aligner.utils.daaMerger.Hit.FrameDirection;
import ella.model.io.MyParameters;

public class HitFilter {
	
	public static ArrayList<Hit> run(ArrayList<Hit> hits, double lambda, double K) {

		ArrayList<Hit> passedHits = new ArrayList<Hit>();
		for (Hit h1 : hits) {
			
			// checking whether h1 is dominated by another hit
			boolean isDominated = false;
			for (Hit h2 : hits) {
				if (!h1.equals(h2)) {
					double coverage = cmpHitCoverage(h1, h2);
					double bitScore1 = cmpBitScore(h1.getRawScore(), lambda, K);
					double bitScore2 = cmpBitScore(h2.getRawScore(), lambda, K);
					if (coverage > MyParameters.MIN_PROPORTION_COVERAGE && MyParameters.MIN_PROPORTION_SCORE * bitScore2 > bitScore1) {
						isDominated = true;
						break;
					}

				}
			}

			// only non-dominated hits are reported
			if (!isDominated)
				passedHits.add(h1);

		}

		return passedHits;

	}

	private static double cmpBitScore(int rawScore, double lambda, double K) {
		return (new Double(rawScore) * lambda - Math.log(K)) / Math.log(2);
	}

	private static double cmpHitCoverage(Hit h1, Hit h2) {

		int[] h1Coord = getQueryCoordinates(h1);
		int[] h2Coord = getQueryCoordinates(h2);

		// checking if overlap exists
		if (h1Coord[0] > h2Coord[1])
			return 0;
		if (h2Coord[0] > h1Coord[1])
			return 0;

		// computing coverage of h1 by h2
		double l = Math.max(h1Coord[0], h2Coord[0]);
		double r = Math.min(h1Coord[1], h2Coord[1]);
		double overlap = r - l + 1.;
		double length1 = h1Coord[1] - h1Coord[0] + 1;
		double coverage = overlap / length1;

		return coverage;
	}

	private static int[] getQueryCoordinates(Hit h) {
		int queryStart = h.getFrame() == FrameDirection.POSITIVE ? h.getQueryStart() : h.getQueryStart() - h.getQueryLength() + 1;
		int queryEnd = queryStart + h.getQueryLength() - 1;
		int[] coord = { queryStart, queryEnd };
		return coord;
	}

}
