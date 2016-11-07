/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.core;

import edu.duke.humangenome.sam.MultiPileup;
import edu.duke.humangenome.sam.IndelPileup;
import edu.duke.humangenome.util.BaseUtils;
import java.util.Map;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.SamLocusIterator;

/**
 *
 * @author Yongzhuang Liu
 */
public class FeatureSelection {

    private final static int FATHER_INDEX = 0;
    private final static int MOTHER_INDEX = 1;
    private final static int OFFSPRING_INDEX = 2;
    private ReferenceSequence referenceSequence;
    private MultiPileup trioPileup;

    public FeatureSelection(ReferenceSequence referenceSequence, MultiPileup trioPileup) {
        this.referenceSequence = referenceSequence;
        this.trioPileup = trioPileup;
    }

    public String extract() {

        String chrom = trioPileup.getReferenceName();
        //System.out.println("chrom"+chrom);
        int pos = trioPileup.getPosition();
        //System.out.println("pos"+pos);
//        char ref = getRefBase();
//        char alt = getMostAltAllele();
        int[] index = new int[]{FATHER_INDEX, MOTHER_INDEX, OFFSPRING_INDEX};

        StringBuilder recordBuilder = new StringBuilder();

        int[] fAlleleCount = new int[2];
        int[] mAlleleCount = new int[2];
        int[] oAlleleCount = new int[2];

        for (int i = 0; i < index.length; i++) {
            
            IndelPileup pileup = new IndelPileup(referenceSequence,trioPileup.getLocusInfo(index[i]),trioPileup.getRefAllele(), trioPileup.getAltAllele());
            int readDepth = pileup.getDepth();

//            int meanRefBaseQuality = pileup.getMeanRefMappingQuality();
//            int meanAltBaseQuality = pileup.getMeanAltMappingQuality();

            int meanRefMappingQuality = pileup.getMeanRefMappingQuality();
            int meanAltMappingQuality = pileup.getMeanAltMappingQuality();

            int meanRefDistanceToThreePrime = pileup.getMeanRefDistanceToThreePrime();
            int meanAltDistanceToThreePrime = pileup.getMeanAltDistanceToThreePrime();

            double refFractionOfMQ0Reads = pileup.getRefFractionOfMQ0Reads();
            double altFractionOfMQ0Reads = pileup.getAltFractionOfMQ0Reads();

            double refFractionOfSoftClippedReads = pileup.getRefFractionOfSoftClippedReads();
            double altFractionOfSoftClippedReads = pileup.getAltFractionOfSoftClippedReads();

            int forwardRef = pileup.getForwardRefCount();
            int reverseRef = pileup.getReverseRefCount();
            int forwardAlt = pileup.getForwardAltCount();
            int reverseAlt = pileup.getReverseAltCount();

            int strandBias = pileup.getStrandBias();

            int refStrandDirection;
            int altStrandDirection;
            if (pileup.getRefStrandDirection()) {
                refStrandDirection = 1;
            } else {
                refStrandDirection = 0;
            }
            if (pileup.getAltStrandDirection()) {
                altStrandDirection = 1;
            } else {
                altStrandDirection = 0;
            }

            if (index[i] == FATHER_INDEX) {
                fAlleleCount[0] = forwardRef + reverseRef;
                fAlleleCount[1] = forwardAlt + reverseAlt;
            }
            if (index[i] == MOTHER_INDEX) {
                mAlleleCount[0] = forwardRef + reverseRef;
                mAlleleCount[1] = forwardAlt + reverseAlt;
            }
            if (index[i] == OFFSPRING_INDEX) {
                oAlleleCount[0] = forwardRef + reverseRef;
                oAlleleCount[1] = forwardAlt + reverseAlt;
            }
            
            double alleleBalance = 0;//·ÀÖ¹³öÏÖ³ýÁã´íÎó
            if( (pileup.getRefAlleleCount() + pileup.getAltAlleleCount()) != 0 )
            	alleleBalance = (double) pileup.getAltAlleleCount() / (pileup.getRefAlleleCount() + pileup.getAltAlleleCount());
            	
            /*
             double meanRefNearbyIndels = (double) pileup.getNearbyIndels(ref) / (forwardRef + reverseRef);
             double meanAltNearbyIndels = (double) pileup.getNearbyIndels(alt) /  (forwardAlt+ reverseAlt);
             double meanRefNearbyMismatches = (double) (pileup.getOverallMismatches(referenceSequence, ref)) / (forwardRef + reverseRef);
             double meanAltNearbyMismatches = (double) (pileup.getOverallMismatches(referenceSequence, alt) - forwardAlt - reverseAlt) / (forwardAlt+ reverseAlt);

             * 
             */


            double meanRefNearbyIndels = pileup.getMeanRefNearbyIndels();
            double meanAltNearbyIndels = pileup.getMeanAltNearbyIndels();
            double meanRefNearbyMismatches = pileup.getMeanRefNearbyMismatches();
            double meanAltNearbyMismatches = pileup.getMeanAltNearbyMismatches();

            recordBuilder.append(alleleBalance + ","  + readDepth + ",");
            recordBuilder.append(meanRefMappingQuality + "," + meanAltMappingQuality + ",");
            recordBuilder.append(meanRefDistanceToThreePrime + "," + meanAltDistanceToThreePrime + ",");
            recordBuilder.append(refFractionOfMQ0Reads + "," + altFractionOfMQ0Reads + ",");
            recordBuilder.append(refFractionOfSoftClippedReads + "," + altFractionOfSoftClippedReads + ",");
            recordBuilder.append(meanRefNearbyMismatches + "," + meanAltNearbyMismatches + ",");
            recordBuilder.append(meanRefNearbyIndels + "," + meanAltNearbyIndels + ",");
            recordBuilder.append(refStrandDirection + "," + altStrandDirection + ",");
            recordBuilder.append(strandBias + ",");
        }
        int pValueOfFatherToOffspring = getPhredPValue(fAlleleCount, oAlleleCount);
        int pValueOfMotherToOffspring = getPhredPValue(mAlleleCount, oAlleleCount);
        recordBuilder.append(pValueOfFatherToOffspring + "," + pValueOfMotherToOffspring);
        return recordBuilder.toString();
    }
    
    
    public int getPhredPValue(int[] pAlleleCount, int[] oAlleleCount) {
        FisherExact fe = new FisherExact(pAlleleCount[0] + pAlleleCount[1] + oAlleleCount[0] + oAlleleCount[1]);
        double p = fe.getCumlativeP(pAlleleCount[0], pAlleleCount[1], oAlleleCount[0], oAlleleCount[1]);
        if (p > 1) {
            p = 1.0;
        }
        return (int) (Math.log10(p) * (-10));
    }
}
