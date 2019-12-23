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
import java.util.concurrent.ConcurrentHashMap;

import ella.model.aligner.utils.SparseString;

public class ReadHits {

	private ConcurrentHashMap<SparseString, ConcurrentHashMap<Integer, ArrayList<Hit>>> hitMap = new ConcurrentHashMap<SparseString, ConcurrentHashMap<Integer, ArrayList<Hit>>>();

	public void add(Hit h, SparseString gi, int frame) {
		if (!hitMap.containsKey(gi))
			hitMap.put(gi, new ConcurrentHashMap<>());
		ConcurrentHashMap<Integer, ArrayList<Hit>> frameMap = hitMap.get(gi);
		if (!frameMap.containsKey(frame))
			frameMap.put(frame, new ArrayList<>());
		ArrayList<Hit> hits = frameMap.get(frame);
		hits.add(h);
	}

	public ConcurrentHashMap<SparseString, ConcurrentHashMap<Integer, ArrayList<Hit>>> getHitMap() {
		return hitMap;
	}

	public ArrayList<Hit> getAllHits() {
		ArrayList<Hit> allHits = new ArrayList<Hit>();
		for (SparseString gi : hitMap.keySet()) {
			for (int frame : hitMap.get(gi).keySet()) {
				for (Hit h : hitMap.get(gi).get(frame)) {
					allHits.add(h);
				}
			}
		}
		return allHits;
	}

	public void freeFrameHits(SparseString gi, int frame) {
		hitMap.get(gi).remove(frame);
	}

	public void freeGiHits(SparseString gi) {
		for (int frame : hitMap.get(gi).keySet())
			freeFrameHits(gi, frame);
		hitMap.remove(gi);
	}
	
	public void print() {
		for (SparseString gi : hitMap.keySet()) {
			System.out.println(">GI: " + gi);
			for (int frame : hitMap.get(gi).keySet()) {
				System.out.println("\t>>Frame: " + frame);
				for (Hit h : hitMap.get(gi).get(frame)) {
					h.print("\t");
				}
			}
		}
	}

}
