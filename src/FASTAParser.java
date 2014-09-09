import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by Nils Henning on 02-09-2014.
 */
public class FASTAParser {

    public static String Parse(String path) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line,result ="";


        while ((line = br.readLine())!= null) {
            if(line.charAt(0)=='>'||line.charAt(0)==';')continue;
            result +=line;
        }

        if(br != null)br.close();
        return result.replaceAll("\\s","").toUpperCase();
    }
}
