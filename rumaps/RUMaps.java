package rumaps;

import java.util.*;

/**
 * This class represents the information that can be attained from the Rutgers University Map.
 * 
 * The RUMaps class is responsible for initializing the network, streets, blocks, and intersections in the map.
 * 
 * You will complete methods to initialize blocks and intersections, calculate block lengths, find reachable intersections,
 * minimize intersections between two points, find the fastest path between two points, and calculate a path's information.
 * 
 * Provided is a Network object that contains all the streets and intersections in the map
 * 
 * @author Vian Miranda
 * @author Anna Lu
 */
public class RUMaps {
    
    private Network rutgers;

    /**
     * **DO NOT MODIFY THIS METHOD**
     * 
     * Constructor for the RUMaps class. Initializes the streets and intersections in the map.
     * For each block in every street, sets the block's length, traffic factor, and traffic value.
     * 
     * @param mapPanel The map panel to display the map
     * @param filename The name of the file containing the street information
     */
    public RUMaps(MapPanel mapPanel, String filename) {
        StdIn.setFile(filename);
        int numIntersections = StdIn.readInt();
        int numStreets = StdIn.readInt();
        StdIn.readLine();
        rutgers = new Network(numIntersections, mapPanel);
        ArrayList<Block> blocks = initializeBlocks(numStreets); 
        initializeIntersections(blocks);

        for (Block block: rutgers.getAdjacencyList()) {
            Block ptr = block;
            while (ptr != null) {
                ptr.setLength(blockLength(ptr));
                ptr.setTrafficFactor(blockTrafficFactor(ptr));
                ptr.setTraffic(blockTraffic(ptr));
                ptr = ptr.getNext();
            }
        }
    }

    /**
     * **DO NOT MODIFY THIS METHOD**
     * 
     * Overloaded constructor for testing.
     * 
     * @param filename The name of the file containing the street information
     */
    public RUMaps(String filename) {
        this(null, filename);
    }

    /**
     * **DO NOT MODIFY THIS METHOD**
     * 
     * Overloaded constructor for testing.
     */
    public RUMaps() { 
        
    }

    /**
     * Initializes all blocks, given a number of streets.
     * the file was opened by the constructor - use StdIn to continue reading the file
     * @param numStreets the number of streets
     * @return an ArrayList of blocks
     */
    public ArrayList<Block> initializeBlocks(int numStreets) {
        // WRITE YOUR CODE HERE
        ArrayList<Block> arrayBlocks = new ArrayList<>(); 
        
        for(int i = 0; i < numStreets; i++){ 
            String streetName = StdIn.readLine(); 
            int numBlocks = StdIn.readInt(); 
            StdIn.readLine(); 

            for(int j = 0; j < numBlocks; j++){ 
                int blockNumber = StdIn.readInt(); 
                StdIn.readLine();
                int numOfPoints = StdIn.readInt(); 
                StdIn.readLine();
                double roadSize = StdIn.readDouble(); 
                StdIn.readLine();
                Block newBlock = new Block(roadSize, streetName, blockNumber); 

                for(int o = 0; o < numOfPoints; o++){ 
                    int x = StdIn.readInt(); 
                    int y = StdIn.readInt(); 
                    StdIn.readLine();
                    Coordinate newCorrdinate = new Coordinate(x, y); 

                    if(o == 0){  
                        newBlock.startPoint(newCorrdinate);
                    }
                    else{
                        newBlock.nextPoint(newCorrdinate);
                    }
                }

                arrayBlocks.add(newBlock); 

            }
        }
        
        return arrayBlocks; // Replace this line, it is provided so the code compiles
    }

    /**
     * This method traverses through each block and finds
     * the block's start and end points to create intersections. 
     * 
     * It then adds intersections as vertices to the "rutgers" graph if
     * they are not already present, and adds UNDIRECTED edges to the adjacency
     * list.
     * 
     * Note that .addEdge(__) ONLY adds edges in one direction (a -> b). 
     */
    public void initializeIntersections(ArrayList<Block> blocks) {
        // WRITE YOUR CODE HERE
        for (Block block : blocks){ 
            ArrayList<Coordinate> coordinateList = block.getCoordinatePoints(); 
            Coordinate first = coordinateList.get(0); 
            Coordinate last = coordinateList.get(coordinateList.size() - 1); 

            int firstIntersection = rutgers.findIntersection(first); 
            int lastIntersection = rutgers.findIntersection(last);            

            if(firstIntersection == -1){ 
                Intersection newObject = new Intersection(first); 
                block.setFirstEndpoint(newObject); 
                rutgers.addIntersection(newObject); 
                firstIntersection = rutgers.findIntersection(first);
            }
            else{ 
                block.setFirstEndpoint(rutgers.getIntersections()[firstIntersection]); 
            }

            if(lastIntersection == -1){ 
                Intersection newerObject = new Intersection(last);
                block.setLastEndpoint(newerObject);
                rutgers.addIntersection(newerObject);
                lastIntersection = rutgers.findIntersection(last);
            }
            else{
                block.setLastEndpoint(rutgers.getIntersections()[lastIntersection]);
            }


            Block blockA = block.copy(); 
            Block blockB = block.copy(); 

            blockB.setFirstEndpoint(rutgers.getIntersections()[lastIntersection]);
            blockB.setLastEndpoint(rutgers.getIntersections()[firstIntersection]);
            
            rutgers.addEdge(firstIntersection, blockA);
            rutgers.addEdge(lastIntersection, blockB);       
        }
     }

    /**
     * Calculates the length of a block by summing the distances between consecutive points for all points in the block.
     * 
     * @param block The block whose length is being calculated
     * @return The total length of the block
     */
    public double blockLength(Block block) {
        // WRITE YOUR CODE HERE 
        ArrayList<Coordinate> points = block.getCoordinatePoints();
        double distance = 0.0;

        for(int i = 0; i < points.size() - 1; i++){
            Coordinate pt1 = points.get(i);
            Coordinate pt2 = points.get(i+1);
            distance = distance + coordinateDistance(pt1, pt2);
        }
      
        return distance; // Replace this line, it is provided so the code compiles
    }

    /**
     * Use a DFS to traverse through blocks, and find the order of intersections
     * traversed starting from a given intersection (as source).
     * 
     * Implement this method recursively, using a helper method.
     */
    public ArrayList<Intersection> reachableIntersections(Intersection source) {
        // WRITE YOUR CODE HERE
        int index = rutgers.findIntersection(source.getCoordinate());

        if(index == -1){
            return new ArrayList<>();
        }

        boolean[] marked = new boolean[rutgers.getIntersections().length];
        ArrayList<Intersection> edgeTo = new ArrayList<>();

        dfsHelperMethod(index, marked, edgeTo);
        return edgeTo; // Replace this line, it is provided so the code compiles
    }

    private void dfsHelperMethod(int s, boolean[] marked, ArrayList<Intersection> edgeTo){
        marked[s] = true; 
        edgeTo.add(rutgers.getIntersections()[s]);

        Block block = rutgers.adj(s);
        
        while(block != null){
            int index1 = rutgers.findIntersection(block.getFirstEndpoint().getCoordinate());
            int index2 = rutgers.findIntersection(block.getLastEndpoint().getCoordinate());

            if(!marked[index1]){
                dfsHelperMethod(index1, marked, edgeTo);
            }

            if(!marked[index2]){
                dfsHelperMethod(index2, marked, edgeTo);
            }
            block = block.getNext();
        }
    }

   
     

    /**
     * Finds and returns the path with the least number of intersections (nodes) from the start to the end intersection.
     * 
     * - If no path exists, return an empty ArrayList.
     * - This graph is large. Find a way to eliminate searching through intersections that have already been visited.
     * 
     * @param start The starting intersection
     * @param end The destination intersection
     * @return The path with the least number of turns, or an empty ArrayList if no path exists
     */
    public ArrayList<Intersection> minimizeIntersections(Intersection start, Intersection end) {
        // WRITE YOUR CODE HERE
        boolean[] marked = new boolean[rutgers.getIntersections().length];
        Intersection[] edgeTo = new Intersection[rutgers.getIntersections().length];
        Queue<Intersection> q = new Queue<>();

        int index = rutgers.findIntersection(start.getCoordinate());
        int index2 = rutgers.findIntersection(end.getCoordinate());

        if(index == -1 || index2 == -1){
            return new ArrayList<>();
        }

        q.enqueue(rutgers.getIntersections()[index]);
        marked[index] = true;
        edgeTo[index] = null;

        while(!q.isEmpty()){
            Intersection intersection = q.dequeue();
            if(intersection.equals(end)){
                break;
            }

            int index3 = rutgers.findIntersection(intersection.getCoordinate());
            Block block = rutgers.adj(index3);

            while(block != null){
                Intersection int1 = block.getFirstEndpoint();
                Intersection int2 = block.getLastEndpoint();

                if(!marked[rutgers.findIntersection(int1.getCoordinate())]){
                    marked[rutgers.findIntersection(int1.getCoordinate())] = true;
                    edgeTo[rutgers.findIntersection(int1.getCoordinate())] = intersection;
                    q.enqueue(int1);
                }

                if(!marked[rutgers.findIntersection(int2.getCoordinate())]){
                    marked[rutgers.findIntersection(int2.getCoordinate())] = true;
                    edgeTo[rutgers.findIntersection(int2.getCoordinate())] = intersection;
                    q.enqueue(int2);
                }

                block = block.getNext();
            }
        }
        ArrayList<Intersection> path = new ArrayList<>();
        Intersection current = end;

        while(current != null){
            path.add(current);
            current = edgeTo[rutgers.findIntersection(current.getCoordinate())];
        }

        Collections.reverse(path);
        return path; // Replace this line, it is provided so the code compiles
    }

    /**
     * Finds the path with the least traffic from the start to the end intersection using a variant of Dijkstra's algorithm.
     * The traffic is calculated as the sum of traffic of the blocks along the path.
     * 
     * What is this variant of Dijkstra?
     * - We are using traffic as a cost - we extract the lowest cost intersection from the fringe.
     * - Once we add the target to the done set, we're done. 
     * 
     * @param start The starting intersection
     * @param end The destination intersection
     * @return The path with the least traffic, or an empty ArrayList if no path exists
     */
    public ArrayList<Intersection> fastestPath(Intersection start, Intersection end) {
        // WRITE YOUR CODE HERE
        
        int index = rutgers.findIntersection(start.getCoordinate());
        int index2 = rutgers.findIntersection(end.getCoordinate());
        if (index == -1 || index2 == -1) {
            return new ArrayList<>();
        }
    
        double[] cost = new double[rutgers.getIntersections().length];
        Arrays.fill(cost, Double.POSITIVE_INFINITY);
        cost[index] = 0;

        PriorityQueue<Intersection> pq = new PriorityQueue<>(Comparator.comparingDouble(i -> cost[rutgers.findIntersection(i.getCoordinate())]));
        pq.add(rutgers.getIntersections()[index]);
    
        Intersection[] edgeTo = new Intersection[rutgers.getIntersections().length];
    
        while (!pq.isEmpty()){
            Intersection intersection = pq.poll();
            int indeX = rutgers.findIntersection(intersection.getCoordinate());
    
            if (intersection.equals(end)){
                break;
            } 
    
            Block block = rutgers.adj(indeX);

            while (block != null){
                Intersection first = block.getFirstEndpoint();
                Intersection second = block.getLastEndpoint();
    
                double firstD = cost[indeX] + block.getTraffic();
                double secondD = cost[indeX] + block.getTraffic();
    
                int integer = rutgers.findIntersection(first.getCoordinate());

                if (firstD < cost[integer]){
                    cost[integer] = firstD;
                    edgeTo[integer] = intersection;
                    pq.add(first);
                }
    
                int integer2 = rutgers.findIntersection(second.getCoordinate());

                if (secondD < cost[integer2]){
                    cost[integer2] = secondD;
                    edgeTo[integer2] = intersection;
                    pq.add(second);
                }
    
                block = block.getNext();
            }
        }
    
        ArrayList<Intersection> path = new ArrayList<>();
        Intersection intersection2 = end;
        while(intersection2 != null){
            path.add(intersection2);
            int indEx = rutgers.findIntersection(intersection2.getCoordinate());
            intersection2 = edgeTo[indEx];
        }
    
        Collections.reverse(path);
        return path;
    }

    /**
     * Calculates the total length, average experienced traffic factor, and total traffic for a given path of blocks.
     * 
     * You're given a list of intersections (vertices); you'll need to find the edge in between each pair.
     * 
     * Compute the average experienced traffic factor by dividing total traffic by total length.
     *  
     * @param path The list of intersections representing the path
     * @return A double array containing the total length, average experienced traffic factor, and total traffic of the path (in that order)
     */
    public double[] pathInformation(ArrayList<Intersection> path) {
        // WRITE YOUR CODE HERE
        double length = 0.0;
        double traffic = 0.0;

        for(int i = 0; i < path.size() - 1; i++){
            Intersection first = path.get(i);
            Intersection last = path.get(i+1);

            int index = rutgers.findIntersection(first.getCoordinate());
            Block block = rutgers.adj(index);

            while(block != null){
                if((block.getFirstEndpoint().equals(first) && block.getLastEndpoint().equals(last)) || (block.getFirstEndpoint().equals(last) && block.getLastEndpoint().equals(first))){
                    length = length + block.getLength();
                    traffic = traffic + block.getTraffic();
                    break;
                }
                block = block.getNext();
            }
        }
        double trafficFactor = 0.0;
        if(length != 0){
            trafficFactor = traffic / length;
        }
        return new double[] {length, trafficFactor, traffic}; // Replace this line, it is provided so the code compiles
    }

    /**
     * Calculates the Euclidean distance between two coordinates.
     * PROVIDED - do not modify
     * 
     * @param a The first coordinate
     * @param b The second coordinate
     * @return The Euclidean distance between the two coordinates
     */
    private double coordinateDistance(Coordinate a, Coordinate b) {
        // PROVIDED METHOD

        double dx = a.getX() - b.getX();
        double dy = a.getY() - b.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }

    /**
     * **DO NOT MODIFY THIS METHOD**
     * 
     * Calculates and returns a randomized traffic factor for the block based on a Gaussian distribution.
     * 
     * This method generates a random traffic factor to simulate varying traffic conditions for each block:
     * - < 1 for good (faster) conditions
     * - = 1 for normal conditions
     * - > 1 for bad (slower) conditions
     * 
     * The traffic factor is generated with a Gaussian distribution centered at 1, with a standard deviation of 0.2.
     * 
     * Constraints:
     * - The traffic factor is capped between a minimum of 0.5 and a maximum of 1.5 to avoid extreme values.
     * 
     * @param block The block for which the traffic factor is calculated
     * @return A randomized traffic factor for the block
     */
    public double blockTrafficFactor(Block block) {
        double rand = StdRandom.gaussian(1, 0.2);
        rand = Math.max(rand, 0.5);
        rand = Math.min(rand, 1.5);
        return rand;
    }

    /**
     * Calculates the traffic on a block by the product of its length and its traffic factor.
     * 
     * @param block The block for which traffic is being calculated
     * @return The calculated traffic value on the block
     */
    public double blockTraffic(Block block) {
        // PROVIDED METHOD
        
        return block.getTrafficFactor() * block.getLength();
    }

    public Network getRutgers() {
        return rutgers;
    }




    
    








}
