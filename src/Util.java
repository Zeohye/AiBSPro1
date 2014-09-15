import java.util.Random;

/**
 * Created by Nils Henning on 02-09-2014.
 */
public class Util {

    public static int max3(int v1, int v2, int v3){
        return Math.max(v1,Math.max(v2,v3));
    }
    public static int min3(int v1, int v2, int v3){
        return Math.min(v1,Math.min(v2,v3));
    }
    public static void printMatrix(int[][] matrix){
        for(int i = 0; i<matrix.length;i++){
            for(int j = 0; j<matrix[0].length;j++){
                System.out.print(matrix[i][j]+"\t");
            }
            System.out.println();
        }
    }

    public static String generateString(int length){
        char[] chars = "ACGT".toCharArray();
        StringBuilder sb = new StringBuilder();

        Random random = new Random();

        for(int i = 0; i < length; i++){
            char c = chars[random.nextInt(chars.length)];
            sb.append(c);
        }


        return sb.toString();
    }
}
