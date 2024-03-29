package ella.model.io.converting;

import ella.model.aligner.utils.ScoringMatrix;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.File;
import java.math.BigInteger;
import java.util.ArrayList;

public class BlastXMLFile {

    private DocumentBuilderFactory documentBuilderFactory;
    private DocumentBuilder documentBuilder;
    private Document document;
    private TransformerFactory transformerFactory;
    private Transformer transformer;
    private Element blastOutput;
    private Element blastOutputIterations;
    private Element blastOutputParam;


    public BlastXMLFile() {
        try {
            documentBuilderFactory = DocumentBuilderFactory.newInstance();
            documentBuilder = documentBuilderFactory.newDocumentBuilder();
            document = documentBuilder.newDocument();

            transformerFactory = TransformerFactory.newInstance();
            transformer = transformerFactory.newTransformer();

            document.setXmlVersion("1.0");
            transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty(OutputKeys.METHOD, "xml");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");

            blastOutput = document.createElement("BlastOutput");
            blastOutput.setAttributeNS("http://www.ncbi.nlm.nih.gov",
                    "xs:schemaLocation", "http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput.xsd");
            blastOutputParam = document.createElement("BlastOutput_param");
            blastOutputIterations = document.createElement("BlastOutput_iterations");

            document.appendChild(blastOutput);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void initializeFile(String program, String version, String reference, String database, String queryId,
                               String queryDef, String queryLen) {
        Element blastOutputProgram = document.createElement("BlastOutput_program");
        blastOutputProgram.setTextContent(program);
        Element blastOutputVersion = document.createElement("BlastOutput_version");
        blastOutputVersion.setTextContent(version);
        Element blastOutputReference = document.createElement("BlastOutput_reference");
        blastOutputReference.setTextContent(reference);
        Element blastOutputDb = document.createElement("BlastOutput_db");
        blastOutputDb.setTextContent(database);
        Element blastOutputQueryId = document.createElement("BlastOutput_query-ID");
        blastOutputQueryId.setTextContent(queryId);
        Element blastOutputQueryDef = document.createElement("BlastOutput_query-def");
        blastOutputQueryDef.setTextContent(queryDef);
        Element blastOutputQueryLen = document.createElement("BlastOutput_query-len");
        blastOutputQueryLen.setTextContent(queryLen);

        blastOutput.appendChild(blastOutputProgram);
        blastOutput.appendChild(blastOutputVersion);
        blastOutput.appendChild(blastOutputReference);
        blastOutput.appendChild(blastOutputDb);
        blastOutput.appendChild(blastOutputQueryId);
        blastOutput.appendChild(blastOutputQueryDef);
        blastOutput.appendChild(blastOutputQueryLen);

        blastOutput.appendChild(blastOutputParam);
        blastOutput.appendChild(blastOutputIterations);
    }

    public void addUsedProgramParameters(String matrix, String expect, String gapOpen, String gapExtend, String filter) {
        Element parameters = document.createElement("Parameters");

        Element parametersMatrix = document.createElement("Parameters_matrix");
        Element parametersExpect = document.createElement("Parameters_expect");
        Element parametersGapOpen = document.createElement("Parameters_gap-open");
        Element parametersGapExtend = document.createElement("Parameters_gap-extend");
        Element parametersFilter = document.createElement("Parameters_filter");

        parametersMatrix.setTextContent(matrix);
        parametersExpect.setTextContent(expect);
        parametersGapOpen.setTextContent(gapOpen);
        parametersGapExtend.setTextContent(gapExtend);
        parametersFilter.setTextContent(filter);

        parameters.appendChild(parametersMatrix);
        parameters.appendChild(parametersExpect);
        parameters.appendChild(parametersGapOpen);
        parameters.appendChild(parametersGapExtend);
        parameters.appendChild(parametersFilter);

        blastOutputParam.appendChild(parameters);
    }

    public void addHitIteration(int iterNum, ArrayList<DaaHit> daaHits, BigInteger dbLetters, long dbSeqs, DaaReader daaReader, DaaConverter converter) {

        Element iteration = document.createElement("Iteration");

        Element iterationIterNum = document.createElement("Iteration_iter-num");
        Element iterationQueryId = document.createElement("Iteration_query-ID");
        Element iterationQueryDef = document.createElement("Iteration_query-def");
        Element iterationQueryLen = document.createElement("Iteration_query-len");

        iterationIterNum.setTextContent(String.valueOf(iterNum));
        iterationQueryId.setTextContent("Query_" + iterNum);
        iterationQueryDef.setTextContent(daaHits.get(0).getQueryName());
        iterationQueryLen.setTextContent(String.valueOf(daaHits.get(0).getTotalQueryLength()));

        Element iterationHits = document.createElement("Iteration_hits");

        int hitIter = 1;
        for (DaaHit daaHit : daaHits) {
            Element hit = document.createElement("Hit");

            Element hitNum = document.createElement("Hit_num");
            Element hitId = document.createElement("Hit_id");
            Element hitDef = document.createElement("Hit_def");
            Element hitAccession = document.createElement("Hit_accession");
            Element hitLen = document.createElement("Hit_len");

            Element hitHsps = document.createElement("Hit_hsps");

            Element hsp = document.createElement("Hsp");

            Element hspNum = document.createElement("Hsp_num");
            Element hspBitScore = document.createElement("Hsp_bit-score");
            Element hspScore = document.createElement("Hsp_score");
            Element hspEvalue = document.createElement("Hsp_evalue");
            Element hspQueryFrom = document.createElement("Hsp_query-from");
            Element hspQueryTo = document.createElement("Hsp_query-to");
            Element hspHitFrom = document.createElement("Hsp_hit-from");
            Element hspHitTo = document.createElement("Hsp_hit-to");
            Element hspQueryFrame = document.createElement("Hsp_query-frame");
            Element hspHitFrame = document.createElement("Hsp_hit-frame");
            Element hspIdentity = document.createElement("Hsp_identity");
            Element hspPositive = document.createElement("Hsp_positive");
            Element hspGaps = document.createElement("Hsp_gaps");
            Element hspAlignLen = document.createElement("Hsp_align-len");
            Element hspQseq = document.createElement("Hsp_qseq");
            Element hspHseq = document.createElement("Hsp_hseq");
            Element hspMidline = document.createElement("Hsp_midline");

            hitNum.setTextContent(String.valueOf(hitIter));
            hitId.setTextContent(daaHit.getReferenceName());
            hitDef.setTextContent("UNKNOWN");
            hitAccession.setTextContent(daaHit.getReferenceName());
            hitLen.setTextContent(String.valueOf(daaHit.getTotalRefLength()));

            hspNum.setTextContent(String.valueOf(1));
            hspBitScore.setTextContent(String.valueOf(converter.getBitScore(daaHit.getRawScore())));
            hspScore.setTextContent(String.valueOf(daaHit.getRawScore()));
            hspEvalue.setTextContent(String.valueOf(converter.getEValue(daaHit.getRawScore(), daaHit.getQueryLength())));
            hspQueryFrom.setTextContent(String.valueOf(daaHit.getQueryStart() + 1));
            if (daaHit.getFrame() < 0) {
                hspQueryTo.setTextContent(String.valueOf(daaHit.getQueryStart() - daaHit.getQueryLength() + 2)); // qend
            } else {
                hspQueryTo.setTextContent(String.valueOf(daaHit.getQueryStart() + daaHit.getQueryLength())); // qend
            }
            hspHitFrom.setTextContent(String.valueOf(daaHit.getRefStart() + 1));
            hspHitTo.setTextContent(String.valueOf(daaHit.getRefStart() + daaHit.getRefLength()));
            hspQueryFrame.setTextContent(String.valueOf(daaHit.getFrame()));
            hspHitFrame.setTextContent(String.valueOf(0));

            String ali1 = daaHit.getAlignment()[0];
            String ali2 = daaHit.getAlignment()[1];
            String matrxType = daaReader.getDaaHeader().getScoreMatrix();
            int gop = daaReader.getDaaHeader().getGapOpen(), gep = daaReader.getDaaHeader().getGapExtend();
            ScoringMatrix scoringMatrix = new ScoringMatrix(matrxType, gop, gep);
            AlignmentStatistics stats = new AlignmentStatistics(ali1, ali2, scoringMatrix);

            hspIdentity.setTextContent(String.valueOf(stats.getNidents()));
            hspPositive.setTextContent(String.valueOf(stats.getNpositives()));

            hspGaps.setTextContent(String.valueOf(stats.getNgaps()));
            hspAlignLen.setTextContent(String.valueOf(ali1.length()));

            hspQseq.setTextContent(ali1);
            hspHseq.setTextContent(ali2);
            hspMidline.setTextContent(String.valueOf(stats.getNidents()));

            hsp.appendChild(hspNum);
            hsp.appendChild(hspBitScore);
            hsp.appendChild(hspScore);
            hsp.appendChild(hspEvalue);
            hsp.appendChild(hspQueryFrom);
            hsp.appendChild(hspQueryTo);
            hsp.appendChild(hspHitFrom);
            hsp.appendChild(hspHitTo);
            hsp.appendChild(hspQueryFrame);
            hsp.appendChild(hspHitFrame);
            hsp.appendChild(hspIdentity);
            hsp.appendChild(hspPositive);
            hsp.appendChild(hspGaps);
            hsp.appendChild(hspAlignLen);
            hsp.appendChild(hspQseq);
            hsp.appendChild(hspHseq);
            hsp.appendChild(hspMidline);

            hitHsps.appendChild(hsp);

            hit.appendChild(hitNum);
            hit.appendChild(hitId);
            hit.appendChild(hitDef);
            hit.appendChild(hitAccession);
            hit.appendChild(hitLen);
            hit.appendChild(hitHsps);

            iterationHits.appendChild(hit);

            hitIter++;
        }

        Element iterationStat = document.createElement("Iteration_stat");

        Element statistics = document.createElement("Statistics");

        Element statisticsDbNum = document.createElement("Statistics_db-num");
        Element statisticsDbLen = document.createElement("Statistics_db-len");
        Element statisticsHspLen = document.createElement("Statistics_hsp-len");
        Element statisticsEffSpace = document.createElement("Statistics_eff-space");
        Element statisticsKappa = document.createElement("Statistics_kappa");
        Element statisticsLambda = document.createElement("Statistics_lambda");
        Element statisticsEntropy = document.createElement("Statistics_entropy");

        statisticsDbNum.setTextContent(String.valueOf(dbSeqs));
        statisticsDbLen.setTextContent(String.valueOf(dbLetters));
        statisticsHspLen.setTextContent(String.valueOf(0));
        statisticsEffSpace.setTextContent(String.valueOf(0));
        statisticsKappa.setTextContent(String.valueOf(converter.K));
        statisticsLambda.setTextContent(String.valueOf(converter.LAMBDA));
        statisticsEntropy.setTextContent(String.valueOf(0));

        statistics.appendChild(statisticsDbNum);
        statistics.appendChild(statisticsDbLen);
        statistics.appendChild(statisticsHspLen);
        statistics.appendChild(statisticsEffSpace);
        statistics.appendChild(statisticsKappa);
        statistics.appendChild(statisticsLambda);
        statistics.appendChild(statisticsEntropy);

        iterationStat.appendChild(statistics);

        iteration.appendChild(iterationIterNum);
        iteration.appendChild(iterationQueryId);
        iteration.appendChild(iterationQueryDef);
        iteration.appendChild(iterationQueryLen);
        iteration.appendChild(iterationHits);
        iteration.appendChild(iterationStat);

        blastOutputIterations.appendChild(iteration);
    }

    public void writeXML(String filepath) throws TransformerException {
        DOMSource domSource = new DOMSource(document);
        StreamResult streamResult = new StreamResult(new File(filepath));
        transformer.transform(domSource, streamResult);
    }
}
