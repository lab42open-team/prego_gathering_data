/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Hashtable;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import static java.lang.Thread.sleep;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;

/**
 *
 * @author Pavlopoulos
 */
public class BringPubMedIDs_for_BioProjects {
   static  Hashtable<String, String> temp_ids = new Hashtable();
    static ArrayList<String> bioprojects = new ArrayList();
    static String codebase = "/data/databases/scripts/gathering_data/bioproject/";
    static String url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=%20bioproject&db=pubmed&id=";
    static StringBuffer buf = new StringBuffer();
static String output_file = "";
    public static void main(String[] args) {
        String input_file = args[0];
        output_file = args[1];
        get_bioprojectIDs(input_file);
        write_to_file(output_file, "");

        for (int i = 0; i < bioprojects.size(); i++) {
            if (i>0 && i % 10 == 0) {
                System.out.println("completed:" + i + "/" + bioprojects.size());
                append_to_file(output_file, buf.toString().replaceAll(",\n", "\n"));
                buf.setLength(0);

            }
            String url_target = url + bioprojects.get(i);
            //url_target = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=%20bioproject&db=pubmed&id=420";
            System.out.println(url_target);
            buf.append(bioprojects.get(i) + "\t");
            try {
                boolean found = false;
                URL url = new URL(url_target);
                URLConnection conn = url.openConnection();

                // open the stream and put it into BufferedReader
                BufferedReader br = new BufferedReader(new InputStreamReader(conn.getInputStream()));

                String inputLine;
                while ((inputLine = br.readLine()) != null) {
                    if (inputLine.contains("pubmed")) {
                        found = true;
                    }
                    if (found == true) {
                        if (inputLine.contains("<Id>")) {
                            String pmid = (inputLine.trim().replaceAll("<Id>", "").replace("</Id>", ""));
                            buf.append(pmid + ",");
                        }
                    }
                }
                br.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
            buf.append("\n");
            
            try {
                sleep(250);
            } catch (Exception e) {
            }
            
           // if(i>0 && i%20==0) break;
        }
        
        //System.out.println(buf.toString().replaceAll(",\n", "\n"));
        append_to_file(output_file, buf.toString().replaceAll(",\n", "\n"));
    }

    static void get_bioprojectIDs(String filename) {
        try {
            FileReader fileReader = new FileReader(filename);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            String line;
		int cnt=0;
            while ((line = bufferedReader.readLine()) != null) {
                if (line.length() > 1) {
				temp_ids.put(line,line);
			
		cnt++;
                    }
            }
            bioprojects = Collections.list(temp_ids.keys());
            System.out.println(bioprojects.size());
            fileReader.close();
        } catch (Exception e) {
        }
    }
    
    static void write_to_file( String filename, String message) {
        BufferedWriter output = null;
        try {
            File file = new File(filename);
            output = new BufferedWriter(new FileWriter(file));
            output.write(message);
            output.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
    static void append_to_file(String fileName, String message) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(fileName, true));
            writer.write(message);
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
