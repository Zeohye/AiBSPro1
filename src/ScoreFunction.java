/**
 * Created by Nils Henning on 02-09-2014.
 */
public interface ScoreFunction {

    public int linearGapCost(int a);
    public int affineGapCost(int a, int b, int k);

}
