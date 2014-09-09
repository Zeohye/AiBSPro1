import java.util.ArrayList;
import java.util.Map;

/**
 * Created by Nils Henning on 09-09-2014.
 */
public class linearCostBackTrack {
    ArrayList<String> Sequences1 = new ArrayList<String>();
    ArrayList<String> Sequences2 = new ArrayList<String>();

    public void backTrack (String seq1, String seq2, int[][] table,Map<String,Integer> matrix,String result1,String result2,boolean full,int k){
        int i = seq1.length();
        int j = seq2.length();
        if(!seq1.equals("") && !seq2.equals("")) {
            //check the subCost in -1,-1
            String value =""+ seq1.charAt(i-1) + seq2.charAt(j-1);
            if (i > 0 && j > 0 && table[i][j] == table[i-1][j-1] + matrix.get(value)) {
                String res1 = seq1.charAt(i-1) + result1;
                String res2 = seq2.charAt(j-1) + result2;
                backTrack(seq1.substring(0, i-1), seq2.substring(0, j-1), table, matrix, res1, res2, full,k);
                if (!full) return;
            }
            //Check Insertion block
            if (table[i][j] == table[i - 1][j] - ScoreFunctions.linierGapCost(k)) {
                String res1=result1;
                String res2=result2;
                res2 = "-" + res2;
                res1 = seq1.charAt(seq1.length()-1) + res1;
                backTrack(seq1.substring(0, i -1), seq2, table, matrix, res1, res2, full,k);
                if (!full) return;
            }
            //Check Deletion block
            if (table[i][j] == table[i][j - 1] - ScoreFunctions.linierGapCost(k)) {
                String res1 = result1;
                String res2 = result2;
                res1 = "-" + res1;
                res2 = seq2.charAt(seq2.length() -1) + res2;
                backTrack(seq1, seq2.substring(0, j - 1), table, matrix, res1, res2, full,k);
                if (!full) return;
            }
            return;
        }
        //reached an empty string or 2.
        //add rest of sequence 1
        if(!seq1.equals("")) {
            for (int n = 0; n < seq2.length(); n++) {
                result1 = "-" + result1;
                result2 = seq2.charAt(seq2.length() - n -1) + result2;
            }
        }
        //add rest of sequence 2
        if(!seq2.equals("")) {
            for (int n = 0; n < seq1.length(); n++) {
                result2 = "-" + result2;
                result1 = seq1.charAt(seq1.length() - n - 1) + result1;
            }
        }
        Sequences1.add(result1);
        Sequences2.add(result2);
    }
    public ArrayList<String> getSequenses1 (){
        return Sequences1;
    }
    public ArrayList<String> getSequenses2 (){
        return Sequences2;
    }

}
