import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by Nils Henning on 02-09-2014.
 */
public class matrixParser {
        public static HashMap Parse(String path) throws IOException {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            HashMap<String,Integer> matrix = new HashMap<String, Integer>();
            String[] header = br.readLine().toUpperCase().split(" ");
            int j=0;
            while ((line = br.readLine())!= null) {
                String[] values = line.split(" ");
                for (int i = 0; i < header.length; i++) {
                        matrix.put(header[i]+header[j],Integer.valueOf(values[i]));
                }
                j++;
            }
            if(br != null)br.close();
            return matrix;
        }
}

