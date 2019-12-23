package ella.model.aligner.aligning.query;

import java.util.ArrayList;
import java.util.HashMap;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.aligning.extending.AlignmentExtender;
import ella.model.aligner.aligning.extending.FrameHit;
import ella.model.aligner.aligning.extending.FrameRun;
import ella.model.aligner.aligning.seeding.BasicSeedFinder;
import ella.model.aligner.aligning.seeding.Seed;
import ella.model.aligner.aligning.seeding.chaining.ChainBuilder;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.indexing.SuffixArrayTable;
import ella.model.aligner.utils.DNAPacker;
import ella.model.aligner.utils.MyStrand;
import ella.model.aligner.utils.SparseString;
import ella.model.aligner.utils.frameshiftAligner.DP_Matrix;
import ella.model.aligner.utils.frameshiftAligner.DP_Matrix_last;
import ella.model.aligner.utils.streams.MySeedStream;
import ella.model.aligner.utils.wrapper.ObjectArrayWrapper;

public class QueryContainer {

    private int rank;
    private ArrayList<FrameRun> frameRuns = new ArrayList<FrameRun>();
    private ArrayList<ReferenceHit> refHits = new ArrayList<ReferenceHit>();
    private HashMap<ObjectArrayWrapper, MySeedStream> plusSeedMap = new HashMap<ObjectArrayWrapper, MySeedStream>();
    private HashMap<ObjectArrayWrapper, MySeedStream> minusSeedMap = new HashMap<ObjectArrayWrapper, MySeedStream>();

    private SparseString id;
    private byte[] packedSequence;
    private int totalQueryLength, bits;

    private Aligner aligner;

    public QueryContainer(SparseString id, byte[] packedSequence, int length, int bits, int rank, Aligner aligner) {
        this.id = id;
        this.packedSequence = packedSequence;
        this.totalQueryLength = length;
        this.bits = bits;
        this.rank = rank;
        this.aligner = aligner;
    }

    public void computeGappedAlignments(IndexText text, DP_Matrix dpMatrix, AlignmentExtender alignmentExtender) {
        String dnaPlus = getStrandSpecificDNA(MyStrand.FWD);
        String dnaMinus = getStrandSpecificDNA(MyStrand.REV);
        for (ReferenceHit refHit : refHits) {
            String query = refHit.getStrand().equals(MyStrand.FWD) ? dnaPlus : dnaMinus;
            ArrayList<FrameRun> runs = new ArrayList<FrameRun>();
            alignmentExtender.run(refHit, text, query, this, dpMatrix, aligner, runs, false);
            frameRuns.addAll(runs);
        }
    }

    public ArrayList<FrameRun> getFrameRuns() {
        return frameRuns;
    }

    public void computeSeedings(IndexText text, SuffixArrayTable saTable, int m, BasicSeedFinder basicSeedFinder) {
        basicSeedFinder.run(this, text, saTable, m, aligner);
        // TopSeedFinder.run(this, text, saTable, m, ella.model.aligner);
    }

    public void addRefHit(ReferenceHit refHit) {
        refHits.add(refHit);
    }

    public void chainSeedings(IndexText text, ChainBuilder chainBuilder) {
        chainBuilder.run(this, text);
        plusSeedMap = null;
        minusSeedMap = null;
    }

    public void addSeed(Seed seed, ObjectArrayWrapper info, MyStrand strand) {
        if (strand == MyStrand.FWD) {
            plusSeedMap.putIfAbsent(info, new MySeedStream());
            plusSeedMap.get(info).add(seed);
        } else {
            minusSeedMap.putIfAbsent(info, new MySeedStream());
            minusSeedMap.get(info).add(seed);
        }
    }

    public String getStrandSpecificDNA(MyStrand strand) {
        char[] sigma = {'A', 'C', 'G', 'T', 'T', 'G', 'C', 'A'};
        StringBuilder buf = new StringBuilder();
        for (byte a : DNAPacker.unpackSequence(packedSequence, totalQueryLength, bits)) {
            int offset = strand.equals(MyStrand.FWD) ? 0 : 4;
            int i = ((int) a) + offset;
            buf.append(sigma[i]);
        }
        if (strand == MyStrand.REV)
            buf.reverse();
        return buf.toString();
    }

    public int getWeight() {
        int weight = 0;
        for (FrameRun run : frameRuns)
            weight += run.getFrameHits().size();
        return weight;
    }


    public HashMap<ObjectArrayWrapper, MySeedStream> getPlusSeedMap() {
        return plusSeedMap;
    }

    public HashMap<ObjectArrayWrapper, MySeedStream> getMinusSeedMap() {
        return minusSeedMap;
    }

    public ArrayList<ReferenceHit> getRefHits() {
        return refHits;
    }

    public byte[] getPackedSequence() {
        return packedSequence;
    }

    public int getTotalQueryLength() {
        return totalQueryLength;
    }

    public SparseString getId() {
        return id;
    }

    public int getRank() {
        return rank;
    }

    public void freeMemory() {
        frameRuns = null;
        refHits = null;
        plusSeedMap = null;
        minusSeedMap = null;
        packedSequence = null;
    }

}
