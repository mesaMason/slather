package slather.g8;

import slather.sim.Cell;
import slather.sim.Point;
import slather.sim.Move;
import slather.sim.Pherome;
import java.util.*;


public class Player implements slather.sim.Player {

    private Random gen;
    private double d;
    private int t;

    private static final int NUMDIRECTIONS = 4;   // CONSTANTS - number of directions
    private static final int DURATION_CAP = 4;   // CONSTANTS - cap on traveling
    private static final int FRIENDLY_AVOID_DIST = 5; // CONSTANTS - for random walker: turn away from friendlies that are in this radius
    private static final double MAX_MOVEMENT = 1; // CONSTANTS - maximum movement rate (assumed to be 1mm?)
    
    public void init(double d, int t) {
	gen = new Random();
	this.d = d;
	this.t = t;
    }

    /*
      Memory byte setup:
      1 bit strategy (currently either square walker or random walker)

      Memory byte for square walker:
      1 bit strategy
      3 bits empty for now
      2 bits duration
      2 bits direction
      
      Memory byte for random walker:
      1 bit strategy
      7 bits previous direction

      Directions:
      0 = North
      1 = East
      2 = South
      3 = West
     */
    public Move play(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {

        Point currentPosition  = player_cell.getPosition();
        int direction = getDirection(memory);
        int duration = getDuration(memory);
        int maxDuration = getMaxDuration(t, NUMDIRECTIONS);
        int strategy = getStrategy(memory);
        Move nextMove = null;
        
	if (player_cell.getDiameter() >= 2) // reproduce whenever possible
	    return new Move(true, (byte) 0b10000000, (byte) 0);
        
        if (strategy == 0) {
            /* Square walker strategy
             */
            duration++;
            if (duration == maxDuration) {
                direction++;
                direction = direction % NUMDIRECTIONS;
                memory = writeDirection(direction, memory);
            }
            memory = writeDuration(duration, memory, maxDuration);

            Point destination = getNewDest(direction);

            /*
            String s = String.format("%8s", Integer.toBinaryString(memory & 0xFF)).replace(' ','0');
            String p = String.format("current x is %f and current y is %f", currentPosition.x, currentPosition.y);
            String des = String.format("dest x is %f and dest y is %f", destination.x, destination.y);
            System.out.println(s);
            System.out.println(p);
            System.out.println(des);
            */
            
            // if collides, try 20 different random directions
            if (collides(player_cell, destination, nearby_cells, nearby_pheromes)) {
                nextMove = null;
                
                for (int i = 0; i < 20; i++) {
                    int arg = gen.nextInt(120);
                    Point vector = extractVectorFromAngle(arg);
                    if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) {
                        nextMove = new Move(vector, memory);
                        break;
                    }
                }
                if (nextMove == null) { // if nothing worked, sit in place
                    nextMove = new Move(new Point(0,0), memory);
                }
            } 
            else { // no collision, free to make square move
                nextMove = new Move(destination, memory);
            }
        } // end square walker strategy
        else {
            /* Random walker strategy:
               Move away from friendly cells that are too close (specified by FRENDLY_AVOID_DIST)
               If no closeby friendly cells to avoid, act like default player (move in straight lines)
             */
            int prevAngle = memory & 0b01111111;

            Iterator<Cell> cell_it = nearby_cells.iterator();
            double sumX = 0;
            double sumY = 0;
            int count = 0;
            Point vector;

            // calculate avg position of nearby friendly cells
            while (cell_it.hasNext()) {
                Cell curr = cell_it.next();
                if (curr.player != player_cell.player)
                    continue;
                if (player_cell.distance(curr) > FRIENDLY_AVOID_DIST) // don't worry about far away friendlies
                    continue;
                Point currPos = curr.getPosition();
                sumX += currPos.x;
                sumY += currPos.y;
                count++;
            }

            if (count == 0) { // case: no friendly cells to move away from
                // if had a previous direction, keep going in that direction
                if (prevAngle > 0) {
                    vector = extractVectorFromAngle( (int)prevAngle);
                    if (!collides( player_cell, vector, nearby_cells, nearby_pheromes)) {
                        nextMove = new Move(vector, memory);
                    }
                }

                // if will collide or didn't have a previous direction, pick a random direction generate move
                if (nextMove == null) {
                    int newAngle = gen.nextInt(120);
                    Point newVector = extractVectorFromAngle(newAngle);
                    byte newMemory = (byte) (0b10000000 | newAngle);
                    nextMove = new Move(newVector, newMemory);
                }
            } else { // case: friendly cells too close, move in opposite direction
                double avgX = sumX / ((double) count);
                double avgY = sumY / ((double) count);

                double towardsAvgX = avgX - currentPosition.x;
                double towardsAvgY = avgY - currentPosition.y;

                double distanceFromAvg = Math.hypot(towardsAvgX, towardsAvgY);
            
                double awayX = (-(towardsAvgX)/distanceFromAvg) * Cell.move_dist;
                double awayY = (-(towardsAvgY)/distanceFromAvg) * Cell.move_dist;

                // clear the previous vector bits
                int newAngle = 0;
                Point newVector = new Point(awayX, awayY);
                byte newMemory = (byte) (0b10000000 | newAngle);
                nextMove = new Move(newVector, newMemory);
            }

            // candidate nextMove written, check for collision
            
            if (collides(player_cell, nextMove.vector, nearby_cells, nearby_pheromes)) {
                nextMove = null;
                int arg = gen.nextInt(120);
                // try 20 times to avoid collision
                for (int i = 0; i < 20; i++) {
                    Point newVector = extractVectorFromAngle(arg);
                    if (!collides(player_cell, newVector, nearby_cells, nearby_pheromes)) {
                        byte newMemory = (byte) (0b10000000 | arg);
                        nextMove = new Move(newVector, newMemory);
                        break;
                    }
                }
                if (nextMove == null) { // if still keeps colliding, stay in place
                    nextMove = new Move(new Point(0,0), (byte) 0b10000000);
                }
            } // end check candidate nextMove collision
        } // end random walker
        
        System.out.println("Next move: " + nextMove.vector.x + ", " + nextMove.vector.y);
        Point estimate = getVector(player_cell, player_cell, nearby_pheromes);
        System.out.println("Estimated last move: " + estimate.x + ", " + estimate.y);
        return nextMove;        
    } // end Move()

    // check if moving player_cell by vector collides with any nearby cell or hostile pherome
    private boolean collides(Cell player_cell, Point vector, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
	Iterator<Cell> cell_it = nearby_cells.iterator();
	Point destination = player_cell.getPosition().move(vector);
	while (cell_it.hasNext()) {
	    Cell other = cell_it.next();
	    if ( destination.distance(other.getPosition()) < 0.5*player_cell.getDiameter() + 0.5*other.getDiameter() + 0.00011)
		return true;
	}
	Iterator<Pherome> pherome_it = nearby_pheromes.iterator();
	while (pherome_it.hasNext()) {
	    Pherome other = pherome_it.next();
            if (other.player != player_cell.player && destination.distance(other.getPosition()) < 0.5*player_cell.getDiameter() + 0.0001)
                return true;
	}
	return false;
    }
    
    // convert an angle (in 3-deg increments) to a vector with magnitude Cell.move_dist (max allowed movement distance)
    private Point extractVectorFromAngle(int arg) {
	double theta = Math.toRadians( 3* (double)arg );
	double dx = Cell.move_dist * Math.cos(theta);
	double dy = Cell.move_dist * Math.sin(theta);
	return new Point(dx, dy);
    }

    private int getDirection(byte mem) {
            return (mem & 3);
    }

    private int getDuration(byte mem) {
            return ((mem >> 2) & 3);
    }
    private byte writeDirection(int direction, byte memory) {
            int actualDirection = direction % NUMDIRECTIONS;
            byte mem = (byte)((memory & 0b11111100) | actualDirection);
            return mem;
    }
    private byte writeDuration(int duration, byte memory, int maxDuration) {
            int actualDuration = duration % maxDuration;
            byte mem = (byte)((memory & 0b11110011) | (actualDuration << 2));
            return mem;
    }

    private int getMaxDuration(int t, int numdirs) {
            return Math.min((t / numdirs), DURATION_CAP);
    }

    private Point getNewDest(int direction) {

            if (direction == 0) {
                    return new Point(0*Cell.move_dist,-1*Cell.move_dist);

            } else if (direction == 1) {
                    return new Point(1*Cell.move_dist,0*Cell.move_dist);

            } else if (direction == 2) {
                    return new Point(0*Cell.move_dist,1*Cell.move_dist);

            } else if (direction == 3) {
                    return new Point(-1*Cell.move_dist,0*Cell.move_dist);

            } else {
                    return new Point(0,0);
            }

    }

    private int getStrategy(byte memory) {
        int strategy = (memory >> 7) & 1;
        return strategy;
    }
    
    /* Estimate the last known direction of a cell given pheromes that can be seen
       Returns as a vector (Point object)
       Method: for given cell, find pherome of the same type that is <= MAX_MOVEMENT
         If more than 1 pherome found, movement cannot be determined, return MAX_MOVEMENT+1,MAX_MOVEMENT+1
         If no pheromes found and cell is at max view distance, movement cannot be determined (otherwise 
           it legitimately did not move!)
     */
    private Point getVector(Cell player_cell, Cell c, Set<Pherome> nearby_pheromes) {
        Iterator<Pherome> pherome_it = nearby_pheromes.iterator();
        double dX;
        double dY;
        int count = 0;
        Pherome closest = null;
        double cRadius = c.getDiameter()/2;
        Point cPos = c.getPosition();
        
        while (pherome_it.hasNext()) {
            Pherome curr = pherome_it.next();
            Point currPos = curr.getPosition();
            if (curr.player != c.player || currPos == cPos)
                continue;
            double distance = c.distance(curr) + cRadius;
            if (distance <= MAX_MOVEMENT) {
                count++;
                closest = curr;
            }
        }

        if (count > 1) { // more than one close pherome, vector cannot be determined
            dX = MAX_MOVEMENT + 1;
            dY = MAX_MOVEMENT + 1;
            System.out.println("Found too many pheromes: " + count);
            return new Point(dX, dY);
        }
        else if (count == 0) { // no pheromes deteced closeby
            double distanceCelltoCell = player_cell.distance(c);
            if ( (distanceCelltoCell + player_cell.getDiameter()/2 + c.getDiameter()/2) >= d ) {
                // other cell is at edge of view, cannot determine vector
                dX = MAX_MOVEMENT + 1;
                dY = MAX_MOVEMENT + 1;
                return new Point(dX, dY);
            } else {
                // other cell well within view and no closeby pheromes, so it actually didn't move
                return new Point(0, 0);
            }
        } 
        else if (count == 1) { // only 1 pherome detected, can get vector
            Point cPosition = c.getPosition();
            Point pPosition = closest.getPosition();
            dX = cPosition.x - pPosition.x;
            dY = cPosition.y - pPosition.y;
            return new Point(dX, dY);
        }
        else {
            dX = MAX_MOVEMENT + 1;
            dY = MAX_MOVEMENT + 1;
            return new Point(dX, dY);
        }
    }
}
