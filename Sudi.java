import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
/**
* PARTNER: CLARA PAKMAN
 */

public class Sudi extends GraphLibrary{

    public static void main(String[] args) throws IOException {
        Graph<Node, Double> sudi = new AdjacencyMapGraph<Node, Double>(); //a graph of nodes (word and count stored together) that has scores as edge labels
        Map<String, Node> hidden = new HashMap<String, Node>(); //a map to navigate the tags and their corresponding words
        BufferedReader inputWords = new BufferedReader(new FileReader("PS5/texts/brown-train-sentences.txt"));
        BufferedReader inputTags = new BufferedReader(new FileReader("PS5/texts/brown-train-tags.txt"));

        System.out.println("> Training the data...");

        train("PS5/texts/brown-train-tags.txt","PS5/texts/brown-train-sentences.txt",sudi,hidden);

        System.out.println("> Testing the sentence: \"my dog is beautiful .\"");

        System.out.println(Viterbi("my dog is beautiful .",sudi,hidden));

        System.out.println("> Percent accuracy for sentence \"my dog is beautiful .\" = ");

        System.out.println(testPercent("my dog is beautiful .", "PRO N V ADJ .", sudi, hidden));

        sudiInteractive(sudi,hidden);

        System.out.println("> Testing the brown test data...");
        System.out.println("* please be patient, this could take a while *");
        String s;
        String sTag;
        Double finalPercent = 0.0;
        while ((s = inputWords.readLine()) != null && (sTag = inputTags.readLine()) != null) {
            finalPercent += testPercent(s,sTag,sudi,hidden);
        }

        System.out.println(finalPercent);

    }

    /**
    * trains a model (observation and transition probabilities) on corresponding lines (sentence and tags) from a pair of training files
    * @param tagsFile -- string of the file containing tags
     * @param wordsFile -- string of the file containing words
     * @param sudi -- graph that will contain all the tags and words and their probabilities
     * @param hidden -- Map of the tags and their corresponding words
     */

    public static void train(String tagsFile, String wordsFile, Graph<Node, Double> sudi, Map<String, Node> hidden) throws IOException{
        String tagLine;
        String wordLine;
        Map<ArrayList<String>, Node> observations = new HashMap<ArrayList<String>, Node>(); // a map of all the words we've seen with their tag that maps to their score
        BufferedReader inputTags = new BufferedReader(new FileReader(tagsFile));
        BufferedReader inputWords = new BufferedReader(new FileReader(wordsFile));

        Node start = new Node("#", 0.0);
        sudi.insertVertex(start);
        hidden.put("#", start);

        while ((tagLine = inputTags.readLine()) != null && (wordLine = inputWords.readLine()) != null) {
            String[] currentTags = tagLine.split(" ");
            String[] currentWords = wordLine.split(" ");

            for (int i = 0; i < currentTags.length; i++) {
                // if we've never seen the tag before, add it to hidden and sudi
                if (!hidden.containsKey(currentTags[i])) {
                    Node curr = new Node(currentTags[i], 0.0);
                    sudi.insertVertex(curr);
                    hidden.put(currentTags[i], curr);
                }

                // store the word and its tag in an arraylist
                ArrayList<String> wordAndTag = new ArrayList<String>();
                wordAndTag.add(currentWords[i]);
                wordAndTag.add(currentTags[i]);

                // if we've never seen that word and tag combo, make it a node in sudi and create an edge to it
                if (!observations.containsKey(wordAndTag)) {
                    Node word = new Node(currentWords[i], 1.0);
                    observations.put(wordAndTag, word);
                    sudi.insertVertex(word);
                    sudi.insertDirected(hidden.get(currentTags[i]), word, 0.0);
                }
                // if we've seen that word and tag combo, increment the in that word's node
                else {
                    observations.get(wordAndTag).count();
                }
                // make a connection from # (start) to the first tag in the sentence in sudi
                if (i == 0) {
                    double count;
                    if (!sudi.hasEdge(hidden.get("#"), hidden.get(currentTags[i]))) {
                        count = 1.0;
                    }
                    else {
                        count = sudi.getLabel(hidden.get("#"), hidden.get(currentTags[i]));
                        count++;
                    }
                    sudi.insertDirected(hidden.get("#"), hidden.get(currentTags[i]), count);
                }
                // make connections between each of the states in the sentence (linearly) in sudi
                else if (i > 0) {
                    double count;
                    if (!sudi.hasEdge(hidden.get(currentTags[i - 1]), hidden.get(currentTags[i]))) {
                        count = 1.0;
                    }
                    else {
                        count = sudi.getLabel(hidden.get(currentTags[i - 1]), hidden.get(currentTags[i]));
                        count++;
                    }
                    sudi.insertDirected(hidden.get(currentTags[i - 1]), hidden.get(currentTags[i]), count);
                }
            }
        }
        // change the the counts to be log scores
        for (String n : hidden.keySet()) {
            int sumLabels = 0;
            int sumWordCounts = 0;
            // find out the sum of all the labels of the outNeighbors of one node
            for (Node outWord : sudi.outNeighbors(hidden.get(n))) {
                if (hidden.get(outWord.word) != null)
                    sumLabels += sudi.getLabel(hidden.get(n), outWord);
                // find out the sum of all the words of one tag
                else
                    sumWordCounts += outWord.count;
            }
            // normalize the counts by their total and log that score
            for (Node outWord : sudi.outNeighbors(hidden.get(n))) {
                if (hidden.get(outWord.word) != null) {
                    double label = sudi.getLabel(hidden.get(n), outWord);
                    sudi.insertDirected(hidden.get(n), outWord, Math.log((label / sumLabels)));
                }
                else {
                    outWord.setCount(Math.log(outWord.count / sumWordCounts));
                }
            }
        }
    }

    /**
     * finds the best sequence of tags for a line (sequence of words)
     * @param s -- string of words to find the tags of
     * @param sudi -- graph that will contain all the tags and words and their probabilities
     * @param hidden -- Map of the tags and their corresponding words
     * @return an arrayList of the tags of the words in the sentence
     */
    public static ArrayList<String> Viterbi(String s, Graph<Node,Double> sudi, Map<String,Node> hidden){
        String[] words = s.split(" ");
        ArrayList<String> currStates = new ArrayList<String>();
        //currStates = { start }
        currStates.add("#");
        Map<String, Double> currScores = new HashMap<String, Double>();
        // currScores = map { start=0 }
        currScores.put("#", 0.0);
        Map<String, ArrayList<String>> backMap = new HashMap<String, ArrayList<String>>();
        // a map of the tag traces from a starting tag
        backMap.put("#", new ArrayList<String>());
        // maxBack = the highest value for one of the states, which is used to pick the right backtrace from the backMap
        String maxBack = null;

        // for i from 0 to # observations - 1
            for (int i = 0; i < words.length; i++) {
                // nextStates = {}
                ArrayList<String> nextStates = new ArrayList<String>();
                // nextScores = empty map
                Map<String, Double> nextScores = new HashMap<String, Double>();
                // the backtrace to be visited next
                Map<String, ArrayList<String>> nextBack = new HashMap<String, ArrayList<String>>();
                for (String currState : currStates) {
                    for (Node nextState : sudi.outNeighbors(hidden.get(currState))) {
                        if (hidden.containsKey(nextState.word)) {
                            if (!nextStates.contains(nextState.word)) {
                                nextStates.add(nextState.word);
                            }

                            double obsScore = -16.0;
                            for (Node word : sudi.outNeighbors(hidden.get(nextState.word))) {
                                if (word.word.equals(words[i])) {
                                    obsScore = word.count;
                                }
                            }
                            if (currScores.get(currState) == null) currScores.put(currState, 0.0);
                            // nextScore = path to here + take a step to there + make the observation there
                            double nextScore = currScores.get(currState) + sudi.getLabel(hidden.get(currState), nextState) + obsScore;

                            if (!nextScores.containsKey(nextState.word) || nextScore > nextScores.get(nextState.word)) {
                                nextScores.put(nextState.word, nextScore);
                                ArrayList<String> back = new ArrayList<String>();
                                back.addAll(backMap.get(currState));
                                back.add(nextState.word);
                                nextBack.put(nextState.word, back); // add to the next backtrace to check
                            }
                        }
                    }
                }

                currStates = nextStates;
                currScores = nextScores;
                backMap = nextBack;
                maxBack = null;
                // finds the tag associated with the highest total score
                double maxBackScore = Double.NEGATIVE_INFINITY;
                for (String str : backMap.keySet()) {
                    if (currScores.get(str) > maxBackScore) {
                        maxBackScore = currScores.get(str);
                        maxBack = str;
                    }
                }
        }
        // returns the highest scoring backtrace
        return backMap.get(maxBack);
    }

    /**
     * allows the user to interact directly with the viterbi method
     * @param sudi -- graph that will contain all the tags and words and their probabilities
     * @param hidden -- Map of the tags and their corresponding words
     * @return an arrayList of the tags of the words in the sentence
     */
    public static void sudiInteractive(Graph<Node,Double> sudi, Map<String,Node> hidden){
        Scanner in = new Scanner(System.in);
        System.out.println("------------------------------");
        System.out.println(" Welcome to the testing mode! ");
        System.out.println("*  please type your sentence *");
        System.out.println("*  * type -1 to exit test *  *");
        System.out.println("-----------------------------");
        System.out.print("> ");
        String s = "";
        while((s=in.nextLine())!=null && (!s.equals("-1"))){
            s = s.toLowerCase();
            System.out.println(Viterbi(s,sudi,hidden));
            System.out.print("> ");
        }
    }

    /**
     * finds the percent accuracy of the training and viterbi algorythms
     * @param s -- string of words to find the tags of
     * @param sTags -- string of tags to check viterbi with
     * @param sudi -- graph that will contain all the tags and words and their probabilities
     * @param hidden -- Map of the tags and their corresponding words
     * @return the percent accuracy between the inputted tags and the tags viterbi returns
     */
    public static double testPercent(String s, String sTags, Graph<Node,Double> sudi, Map<String,Node> hidden){
        double rightWords = 0;
        double wrongWords = 0;
        String[] tags = sTags.split(" ");
        ArrayList<String> outTags = Viterbi(s,sudi,hidden);
        for(int c = 0; c < outTags.size(); c++){
            if(outTags.get(c).equals(tags[c])) rightWords++;
            else wrongWords++;
        }
        return (rightWords/(rightWords+wrongWords))*100;
    }

}