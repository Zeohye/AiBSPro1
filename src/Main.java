import java.io.IOException;
import java.util.Map;

/**
 * Created by Nils Henning on 02-09-2014.
 */
public class Main {

    static int[][] scoreTable;
    static int[][] insTable;
    static int[][] delTable;
    public static void main(String[] args){
        String seq1,seq2;
        int a = Integer.parseInt(args[0]);
        int b = Integer.parseInt(args[1]);
        Map matrix = null;
        try {
            seq1 = FASTAParser.Parse("input/seq1.fasta");
        } catch (IOException e) {
            seq1 = "";
        }
        try {
            seq2 = FASTAParser.Parse("input/seq2.fasta");
        } catch (IOException e) {
            seq2="";
        }

        try {
            matrix = matrixParser.Parse("input/scoreMatrix.txt");
        } catch (IOException e) {
            e.printStackTrace();
        }
       // Util.printMatrix(fillTableLinear(seq1, seq2, matrix, 5));
        fillTableLinear(seq1, seq2, matrix, 5);
        linearCostBackTrack lCallBack = new linearCostBackTrack();
        lCallBack.backTrack(seq1,seq2,scoreTable,matrix,"","",true,a);

        System.out.println("Linear solutions");
        System.out.println(">seq1");
        System.out.println(lCallBack.getSequenses1().get(0));
        System.out.println(">seq2");
        System.out.println(lCallBack.getSequenses2().get(0));
        System.out.println(">Number of optimal alignments");
        System.out.println(lCallBack.getSequenses1().size());
        System.out.println(">Optimal Score");
        System.out.println(scoreTable[seq1.length()][seq2.length()]);

        //Util.printMatrix(fillTableAffine(seq1, seq2, matrix, a, b));
        fillTableAffine(seq1, seq2, matrix, a, b);
        affineCostBackTrack callback = new affineCostBackTrack();
        callback.backTrack(seq1,seq2,scoreTable,matrix,"","",true,a,b);

        System.out.println("Affine solutions");
        System.out.println(">seq1");
        System.out.println(callback.getSequenses1().get(0));
        System.out.println(">seq2");
        System.out.println(callback.getSequenses2().get(0));
        System.out.println(">Number of optimal alignments");
        System.out.println(callback.getSequenses1().size());
        System.out.println(">Optimal Score");
        System.out.println(scoreTable[seq1.length()][seq2.length()]);
    }

    public static int[][] fillTableLinear(String seq1, String seq2, Map<String, Integer> matrix, int gapCost){
        scoreTable = new int [seq1.length()+1][seq2.length()+1];
        scoreTable[0][0]=0;
        for(int i = 1; i<=seq1.length();i++){
            scoreTable [i][0]=-1*i*gapCost;
        }
        for(int j = 1; j<=seq2.length(); j++){
            scoreTable [0][j]=-1*j*gapCost;
        }

        for(int i = 1; i<=seq1.length();i++){
            for(int j = 1; j<=seq2.length();j++){
                int v1,v2,v3;
                String value = ""+seq1.charAt(i-1)+seq2.charAt(j-1);
                v1 = scoreTable[i-1][j-1]+matrix.get(value);
                v2 = scoreTable[i-1][j]-gapCost;
                v3 = scoreTable[i][j-1]-gapCost;
                scoreTable[i][j] = Util.max3(v1, v2, v3);
            }
        }
        return scoreTable;
    }

    public static int[][] fillTableAffine(String seq1, String seq2, Map<String, Integer> matrix, int a, int b){

        scoreTable = new int[seq1.length()+1][seq2.length()+1];
        delTable = new int[seq1.length()+1][seq2.length()+1];
        insTable = new int[seq1.length()+1][seq2.length()+1];
        for(int i =0;i<=seq1.length();i++){
            for(int j =0; j<=seq2.length();j++){
                if(i==0&&j==0) {
                    scoreTable[0][0] = 0;
                    continue;
                }
                //delTable
                if(i==1)
                    delTable[i][j]=scoreTable[i - 1][j] - (a + b);
                else if(i>1)
                    delTable[i][j]=Math.max(scoreTable[i - 1][j] - (a + b), delTable[i - 1][j] - a);
                //insTable
                if(j==1)
                    insTable[i][j]=scoreTable[i][j-1]-(a+b);
                else if(j>1)
                    insTable[i][j]=Math.max(scoreTable[i][j-1]-(a+b),insTable[i][j-1]-a);

                int v1,v2,v3;
                if(i>0 && j>0) {
                    String value = "" + seq1.charAt(i - 1) + seq2.charAt(j - 1);
                    v1 = scoreTable[i - 1][j - 1] + (matrix.get(value));
                }else
                    v1 = Integer.MIN_VALUE;

                if(i>0)
                    v2 = delTable[i][j];
                else
                    v2 = Integer.MIN_VALUE;

                if(j>0)
                    v3 = insTable[i][j];
                else
                    v3 = Integer.MIN_VALUE;

                scoreTable[i][j] = Util.max3(v1, v2, v3);
            }
        }
        return scoreTable;
    }
}
