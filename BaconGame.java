import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class BaconGame extends GraphLibrary{

    public static void main(String[] args) throws FileNotFoundException, IOException {
        Map<Integer,String> actors = new HashMap<Integer, String>();
        Map<Integer,String> movies = new HashMap<Integer, String>();
        Map<Integer,ArrayList<Integer>> moviesToActors = new HashMap<Integer, ArrayList<Integer>>();
        Graph<String, Set<String>> g = new AdjacencyMapGraph<String, Set<String>>();
        GraphLibrary gLib = new GraphLibrary();

        /*
        * reading in the actors and their ID's using an array to differentiate between the ID, the | symbol, and the actor name
         */
        String filename = "actorsTest";
        BufferedReader input = new BufferedReader(new FileReader("PS4/"+filename + ".txt"));
        String line = "";
        while((line = input.readLine()) != null){
            String[] elements = line.split("\\|");
//            elements[0] = line.substring(0,line.indexOf("|"));
//            elements[1] = line.substring(line.indexOf("|")+1);
            actors.put(Integer.parseInt(elements[0]),elements[1]);
            g.insertVertex(elements[1]);
        }

        /*
         * reading in the movies and their ID's using an array to differentiate between the ID, the | symbol, and the movie name
         */
        filename = "moviesTest";
        input = new BufferedReader(new FileReader("PS4/"+filename + ".txt"));
        line = "";
        while((line = input.readLine()) != null){
            String[] elements = line.split("\\|");
//            elements[0] = line.substring(0,line.indexOf("|"));
//            elements[1] = line.substring(line.indexOf("|")+1);
            movies.put(Integer.parseInt(elements[0]),elements[1]);
        }

        /*
         * reading in the movie ID's and the actor ID's using an array to differentiate between the ID's and the | symbol
         * using a list to keep track of when actors costar in multiple movies
         */
        filename = "movie-actorsTest";
        input = new BufferedReader(new FileReader("PS4/"+filename + ".txt"));
        line = "";
        Set<Integer> moviesSeen = new HashSet<Integer>();
        ArrayList<Integer> actorsInMovie = new ArrayList<Integer>();
        while((line = input.readLine()) != null){
            String[] elements = line.split("\\|");
//            elements[0] = line.substring(0,line.indexOf("|"));
//            elements[1] = line.substring(line.indexOf("|")+1);
            int movie = Integer.parseInt(elements[0]);
            int actor = Integer.parseInt(elements[1]);
            if(!moviesSeen.contains(movie)){
                actorsInMovie = new ArrayList<Integer>();
                actorsInMovie.add(actor);
                moviesToActors.put(movie,actorsInMovie);
                moviesSeen.add(movie);
            }
            else{
                actorsInMovie.add(actor);
            }
        }

        /*
        * creating a graph of all the actors with all their costarring movies as labels
        *
        * The set called interim collects all the movies that actors with existing edges have costarred in that are currently marked on the graph
         */
        for(int i: moviesToActors.keySet()){
           for(int j = 0; j < moviesToActors.get(i).size(); j++){
               for(int p = 0; p <  moviesToActors.get(i).size(); p++){
                   if(p != j) {
                       String actor1 = actors.get(moviesToActors.get(i).get(j));
                       String actor2 = actors.get(moviesToActors.get(i).get(p));
                       Set<String> costarred;
                       if(g.hasEdge(actor1, actor2)){
                           Set interim = g.getLabel(actor1, actor2);
                           interim.add(movies.get(i));
                           g.insertUndirected(actor1, actor2, interim);
                       }
                       else {
                           costarred = new HashSet<String>();
                           costarred.add(movies.get(i));
                           g.insertUndirected(actor1, actor2, costarred);
                       }
                   }
               }
           }
        }

        /*
        * THE GAME (interface) BEGINS
         */

        String center = "Kevin Bacon";
        boolean quit = false;

        System.out.println("Welcome to the Kevin Bacon Game! \nControls: \nc <#>: list top (positive number) or bottom (negative) <#> centers of the universe, sorted by average separation\n" +
                "d <low> <high>: list actors sorted by degree, with degree between low and high\n" +
                "i: list actors with infinite separation from the current center\n" +
                "p <name>: find path from <name> to current center of the universe\n" +
                "s <low> <high>: list actors sorted by non-infinite separation from the current center, with separation between low and high\n" +
                "u <name>: make <name> the center of the universe\n" +
                "q: quit game");

        Scanner in = new Scanner(System.in);
        while(quit!=true) {
            String s = in.nextLine();
            /*
            * the input of c allows the user to choose a low bound or a high bound which determines how many actors (least connected or most connected) they'd like to return
            *
            * here we loop through all the actors using the center as an initial grounding root for the graph but then
            * we loop through all the actors as centers for the graph and measure and record their levels of average separation
            * we sort the list in ascneding or descending order depending upon the sign of the user input
            * we then copy all top values (based on the user input) into a list that we print out
            * */
            if (s.equals("c")) {
                System.out.println("Input a positive or negative number depending upon whether you'd like the top most connected actors or the bottom most connected actors");
                int c = in.nextInt();
                Map<String,Double> avSepSort = new HashMap<String,Double>();
                ArrayList<String> allActors = new ArrayList<String>();
                ArrayList<String> topActors = new ArrayList<String>();
                for(String actor: gLib.BFS(g,center).vertices()){
                    allActors.add(actor);
                    for(String actor2: gLib.BFS(g,actor).vertices()) {
                        avSepSort.put(actor, gLib.averageSeparation(gLib.BFS(g, actor), actor2));
                    }
                }
                if(c<0) allActors.sort((String s1, String s2)-> (int)(avSepSort.get(s2)-avSepSort.get(s1)));
                else allActors.sort((String s1, String s2)-> (int)(avSepSort.get(s1)-avSepSort.get(s2)));
                for(int i = 0; i<Math.abs(c); i++){
                    topActors.add(allActors.get(i));
                }
                System.out.println(topActors);
            }

            /*
            * here the user can get a sorted list of degrees of separation from the center
            *
            * we use the upper and lower bound inputs to add specific actors to a list,
            * sorting them by using their inDegrees within the overall graph
             */
            else if (s.equals("d")) {
                ArrayList<Integer> d = new ArrayList<Integer>();
                System.out.println("Input your lower bound");
                d.add(in.nextInt());
                System.out.println("Input your higher bound");
                d.add(in.nextInt());
                System.out.println("Using " + center + " as the center of the universe");
                ArrayList<String> degreeSort = new ArrayList<String>();
                for (String v : verticesByInDegree(g)) {
                    if (g.inDegree(v) >= d.get(0) && g.inDegree(v) <= d.get(1)) {
                        degreeSort.add(v);
                    }
                }
                System.out.println(degreeSort);
            }

            /*
            * here the user can request an output of all the missing vertices (vertices that aren't connected to the center)
            *
            * we do this using our missingVertices method (see GraphLibrary)
             */
            else if (s.equals("i")) {
                System.out.println("Listing the actors with infinite separation from "+center);
                System.out.println(gLib.missingVertices(g, gLib.BFS(g, center)));
            }

            /*
            * here the user can input the name of an actor they would like to trace to the center
            *
            * we do this using our getPath method (see GraphLibrary) and some string manipulation
             */
            else if (s.equals("p")) {
                    System.out.println("Input the name of the actor you'd like to trace to the center");
                    String name = in.nextLine();
                try {
                    if (!gLib.missingVertices(g, gLib.BFS(g, center)).contains(name)) {
                        ArrayList<String> path = (ArrayList<String>) gLib.getPath(gLib.BFS(g, center), name);
                        System.out.println(name + "'s number is " + (path.size() - 1));
                        for (int i = 0; i < path.size() - 1; i++) {
                            System.out.println(path.get(i) + " was in " + g.getLabel(path.get(i), path.get(i + 1)) + " with " + path.get(i + 1));
                        }
                    } else System.out.println("That actor is not connected to " + center);
                }catch(Exception e){
                   System.out.println(name+" does not exist in the universe");
                }
            }

            /*
            * here the user can request an output of sorted actors by separation within an upper and lower bounds
            *
            * we do this using the getPath method (see GraphLibrary) and some comparators with anonymous functions
             */
            else if (s.equals("s")) {
                ArrayList<Integer> sList = new ArrayList<Integer>();
                System.out.println("Input your lower bound");
                sList.add(in.nextInt());
                System.out.println("Input your higher bound");
                sList.add(in.nextInt());
                System.out.println("Using " + center + " as the center of the universe");
                Graph<String, Set<String>> bfs = gLib.BFS(g, center);
                ArrayList<String> degreeSort = new ArrayList<String>();
                for (String v : gLib.BFS(g, center).vertices()) {
                    if (gLib.getPath(bfs, v).size() - 1 >= sList.get(0) && gLib.getPath(bfs, v).size() - 1 <= sList.get(1)) {
                        degreeSort.add(v);
                    }
                }
                degreeSort.sort((String s1, String s2) -> gLib.getPath(bfs, s2).size() - gLib.getPath(bfs, s1).size());
                System.out.println(degreeSort);
            }

            /*
            * here the user can change the center of the universe
             */
            else if (s.equals("u")) {
                String userInput = in.nextLine();
                if(actors.containsValue(userInput)){
                    System.out.println(userInput+" is now the center of the universe!");
                    center = userInput;
                }
                else System.out.println("That actor is not in the universe");
            }

            /*
            * here the user can quit the game
             */
            else if(s.equals("q")){
                System.out.println("Thanks for playing the Kevin Bacon Game!");
                quit = true;
            }
        }

    }



}
