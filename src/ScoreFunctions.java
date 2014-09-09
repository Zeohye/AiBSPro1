/**
 * Created by Nils Henning on 02-09-2014.
 */
public  class ScoreFunctions {

    public static int linierGapCost(int k){
        return -1*k;
    }
    public static int affineGapCost(int a, int b, int k){
        return a*k+b;
    }
}
