package ella.model.io;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class SubjectManager {

	private ConcurrentHashMap<Integer, Integer> subjectMap = new ConcurrentHashMap<Integer, Integer>();
	private AtomicInteger counter = new AtomicInteger(0);

	public int addSubject(int textIndex) {
		if (!subjectMap.containsKey(textIndex))
			subjectMap.put(textIndex, counter.getAndIncrement());
		return subjectMap.get(textIndex);
	}

	public int getSubjectOrder(int textIndex) {
		return subjectMap.get(textIndex);
	}

	public ArrayList<MySubject> getOrderedSubjectList() {
		ArrayList<MySubject> subjectList = new ArrayList<MySubject>();
		for (int key : subjectMap.keySet())
			subjectList.add(new MySubject(key, subjectMap.get(key)));
		Collections.sort(subjectList, new SubjectComparator());
		return subjectList;
	}

	public class MySubject {

		private int locationIndex, order;

		public MySubject(int textIndex, int order) {
			super();
			this.locationIndex = textIndex;
			this.order = order;
		}

		public boolean equals(Object o) {
			if (o instanceof MySubject)
				return locationIndex == ((MySubject) o).getLocationIndex();
			return false;
		}

		public int getLocationIndex() {
			return locationIndex;
		}

		public int getOrder() {
			return order;
		}

		public void setOrder(int order) {
			this.order = order;
		}

	}

	public class SubjectComparator implements Comparator<MySubject> {

		@Override
		public int compare(MySubject s1, MySubject s2) {
			return Integer.compare(s1.getOrder(), s2.getOrder());
		}

	}

}
