#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <inttypes.h>
#include <x86intrin.h>
#include <locale.h>


//  x is after, y is before
int timespec_subtract (struct timespec *result, struct timespec *x, struct timespec *y)
{

    //printf("result2: %lld.%09ldsec taken\n", result2.tv_sec, result2.tv_nsec);

    /* Perform the carry for the later subtraction by updating y. */
    if (x->tv_nsec < y->tv_nsec) {
      int nsec = (y->tv_nsec - x->tv_nsec) / 1000000000 + 1;
      y->tv_nsec -= 1000000000 * nsec;
      y->tv_sec += nsec;
    }
    if (x->tv_nsec - y->tv_nsec > 1000000000) {
      int nsec = (x->tv_nsec - y->tv_nsec) / 1000000000;
      y->tv_nsec += 1000000000 * nsec;
      y->tv_sec -= nsec;
    }

    /* Compute the time remaining to wait.
       tv_usec is certainly positive. */
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_nsec = x->tv_nsec - y->tv_nsec;

    /* Return 1 if result is negative. */
    return x->tv_sec < y->tv_sec;
}

void timespec_add (struct timespec *result, struct timespec *x)
{

    //printf("result2: %lld.%09ldsec taken\n", result2.tv_sec, result2.tv_nsec);

    result->tv_sec = x->tv_sec + result->tv_sec + ((x->tv_nsec + result->tv_nsec) / 1000000000);
    result->tv_nsec = (x->tv_nsec + result->tv_nsec) % 1000000000;

}

void timespec_min_special(struct timespec* timediff, struct timespec* min, int i, int j) {
    if(i == 0 && j == 0) {
        min->tv_sec = timediff->tv_sec;
        min->tv_nsec = timediff->tv_nsec;
    }
    else {
        if(timediff->tv_sec == min->tv_sec) {
            if(timediff->tv_nsec < min->tv_nsec) {
                min->tv_sec = timediff->tv_sec;
                min->tv_nsec = timediff->tv_nsec;
            }
        }
        else if(timediff->tv_sec < min->tv_sec) {
            min->tv_sec = timediff->tv_sec;
            min->tv_nsec = timediff->tv_nsec;
        }
    }
}

void timespec_max_special(struct timespec* timediff, struct timespec* max, int i, int j) {
    if(i == 0 && j == 0) {
        max->tv_sec = timediff->tv_sec;
        max->tv_nsec = timediff->tv_nsec;
    }
    else {
        if(timediff->tv_sec == max->tv_sec) {
            if(timediff->tv_nsec > max->tv_nsec) {
                max->tv_sec = timediff->tv_sec;
                max->tv_nsec = timediff->tv_nsec;
            }
        }
        else if(timediff->tv_sec > max->tv_sec) {
            max->tv_sec = timediff->tv_sec;
            max->tv_nsec = timediff->tv_nsec;
        }
    }
}

uint64_t zxorshiftstate = 0;

uint8_t** grid;
uint32_t* width;
uint32_t* height;
uint32_t total_amount_of_levels = 0;
uint8_t alreadyAllocated = 0;  


void setseed(uint64_t num)
{
    zxorshiftstate = num;
}

uint32_t random()
{
    uint64_t x = zxorshiftstate;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    uint32_t result = (uint32_t)((x * (uint64_t)0x2545F4914F6CDD1DULL) >> 32);
    zxorshiftstate = x;
    return result;
}

typedef struct {
    float pos_x;
    float pos_y;
    float dir_x;
    float dir_y;
} ray_info;

typedef enum {
    NON_OCCUPIED = (uint8_t)0,
    OCCUPIED = (uint8_t)1,
    OUTSIDE = (uint8_t)2,
    PREMATURE_EXIT = (uint8_t)3
} GridResult;

uint8_t readGrid(int32_t x, int32_t y, int32_t level) {
    if(x >= width[level] || x < 0 || 
       y >= height[level] || y < 0) {
       
        return OUTSIDE;
    }

    if(grid[level][x + width[level] * y] == 1) {
        return OCCUPIED;
    }
    else {
        return NON_OCCUPIED;
    }
}

void makeNewOctreeLevel(uint32_t currentLevel) {
    uint32_t previousLevel = currentLevel-1u;
    for (uint32_t i = 0; i < height[currentLevel]; i++)
    {
        for (uint32_t j = 0; j < width[currentLevel]; j++)
        {
            int32_t ll[2] = {j*2,i*2}; //lower left
            int32_t lr[2] = {(j*2)+1,i*2};  //lower right
            int32_t ul[2] = {j*2,(i*2)+1};  //upper left
            int32_t ur[2] = {(j*2)+1,(i*2)+1};  //upper right

            uint8_t occupancy = readGrid(ll[0],ll[1],previousLevel) == OCCUPIED;
            occupancy = (readGrid(lr[0],lr[1],previousLevel) == OCCUPIED) | occupancy;
            occupancy = (readGrid(ul[0],ul[1],previousLevel) == OCCUPIED) | occupancy;
            occupancy = (readGrid(ur[0],ur[1],previousLevel) == OCCUPIED) | occupancy;

            grid[currentLevel][j + i * width[currentLevel]] = occupancy;

        }
    }
    
}

int32_t readOctree(float cellPosition[2]) {
    int32_t cellPositionInt[2] = {(int32_t)floorf(cellPosition[0]),(int32_t)floorf(cellPosition[1])};

    int32_t maxLevel = -1;
    for (int32_t i = 0; i < total_amount_of_levels; i++)
    {
        if(readGrid(cellPositionInt[0] >> i, cellPositionInt[1] >> i,i) == 1) {
            break;
        }
        maxLevel++;
    }

    return maxLevel;
}

typedef struct
{
    uint8_t hit;
    float endposition[2];
    float normal[2];
    float t;
} TraversalOutput;

void printOutput(TraversalOutput output, char methodName[]) {
    printf("%s output:\n", methodName);

    printf("end %.0f %.0f\n", output.endposition[0], output.endposition[1]);
    printf("dir %.0f %.0f\n", output.normal[0], output.normal[1]);
    printf("t: %f\n", output.t);
    if(output.hit) {
        printf("hit true\n");
    }
    else {
        printf("hit false\n");
    }
    
}

typedef struct 
{
    uint64_t correctPositionChecks;
    uint64_t amountPositionChecks;
    uint64_t correctNormalChecks;
    uint64_t amountNormalChecks;

} MethodScore;

typedef struct
{
    float rayPosition[2];
    float rayDir[2];
    float invDir[2];
    float absDir[2];
    float isDirPos[2];
    float isDirNeg[2];
    float isDirNegAlt[2];
    float signDir[2];
} RayInput;

typedef struct
{
    int hitOrNot;
    float cellPosition[2];
    float enteringDirection[2];
    float exitingDirection[2];
    float startT;
    float endT;
} AABBIntersectOutput;

typedef struct 
{
    float axisT;
    uint8_t whichAxis;
    uint8_t ID;
} t_value_axis;

uint8_t compareTraversalOutput(TraversalOutput DDAoutput, TraversalOutput HOutput, MethodScore* HMethod) {
    uint8_t failure = 0;
    uint8_t checkNormal = 0;
    HMethod->amountPositionChecks++;
    if(DDAoutput.hit == 0 && HOutput.hit== 0) {
        HMethod->correctPositionChecks++;
        HMethod->amountNormalChecks++;
        HMethod->correctNormalChecks++;
    }
    else {
        if(DDAoutput.hit == HOutput.hit && DDAoutput.endposition[0] == HOutput.endposition[0] && DDAoutput.endposition[1] == HOutput.endposition[1]) {
            HMethod->correctPositionChecks++;
            checkNormal = 1;
        }
        else {
            failure |= 0b1;
        }

        if(checkNormal) {
            HMethod->amountNormalChecks++;
            if(DDAoutput.normal[0] == HOutput.normal[0] && DDAoutput.normal[1] == HOutput.normal[1]) {
            HMethod->correctNormalChecks++;
            }
            else {
                failure |= 0b10;
            }
        }
    }
    

   

    return failure;
}


void print_t_values(t_value_axis t_values[5]) {
    printf("t_values:\n");
    for (uint8_t i = 0; i < 5; i++)
    {
        printf("ID: %u P: %u T: %f\n", t_values[i].ID,t_values[i].whichAxis,t_values[i].axisT);
    }
    
}


TraversalOutput hero4Total(RayInput ray, AABBIntersectOutput startPosition) {
    TraversalOutput output = {0};

    float cellPosition[2] = {0.0f,0.0f};
    t_value_axis t;
    t.whichAxis = (uint8_t)startPosition.enteringDirection[1];
    t.axisT = startPosition.startT;

    uint64_t closestPowerofTwoCeil[2] = {0};

    for (uint64_t i = 0; i < 32; i++)
    {
        if(width[0] <= (1 << i)) {
            closestPowerofTwoCeil[0] = (1 << i);
            break;
        }
    }
    for (uint64_t i = 0; i < 32; i++)
    {
        if(height[0] <= (1 << i)) {
            closestPowerofTwoCeil[1] = (1 << i);
            break;
        }
    }

    float start_upper_bounds[2] = {closestPowerofTwoCeil[0],closestPowerofTwoCeil[0]};

    float current_upper_bounds[2] = {start_upper_bounds[0], start_upper_bounds[1]};
    float current_lower_bounds[2] = {0};

    int32_t current_level = total_amount_of_levels-1;

    while(!(t.axisT > startPosition.endT || (t.axisT == startPosition.endT && t.whichAxis == (uint8_t)startPosition.exitingDirection[1]))) 
    {
        t_value_axis t_values[5] = {0};

        // x middle plane
        float xplane = (current_lower_bounds[0] + current_upper_bounds[0])/2.0f;
        t_values[0].axisT = (xplane - ray.rayPosition[0]) * ray.invDir[0];
        t_values[0].whichAxis = 0;
        t_values[0].ID = 0;

        // y middle plane
        float yplane = (current_lower_bounds[1] + current_upper_bounds[1])/2.0f;
        t_values[1].axisT = (yplane - ray.rayPosition[1]) * ray.invDir[1];
        t_values[1].whichAxis = 1;
        t_values[1].ID = 1;

        // current t 
        t_values[2].axisT = t.axisT;
        t_values[2].whichAxis = t.whichAxis;
        t_values[2].ID = 2;


        // backside of x plane box
        float xexitplane = ray.isDirNeg[0] ? current_lower_bounds[0] : current_upper_bounds[0];

        t_values[3].axisT = (xexitplane - ray.rayPosition[0]) * ray.invDir[0];
        t_values[3].whichAxis = 0;
        t_values[3].ID = 3;


        // backside of y plane box
        float yexitplane = ray.isDirNeg[1] ? current_lower_bounds[1] : current_upper_bounds[1]; 

        t_values[4].axisT = (yexitplane - ray.rayPosition[1]) * ray.invDir[1];
        t_values[4].whichAxis = 1;
        t_values[4].ID = 3;       


        float inital_position[2] = {
            (xplane - ray.rayPosition[0]) <= 0.0f ? 1.0f : 0.0,
            (yplane - ray.rayPosition[1]) <= 0.0f ? 1.0f : 0.0
        };

        if((xplane - ray.rayPosition[0]) == 0.0f) {
            inital_position[0] = ray.isDirNegAlt[0];
        }
        if((yplane - ray.rayPosition[1]) == 0.0f) {
            inital_position[1] = ray.isDirNegAlt[1];
        }

        for(int16_t i = 2; i < 5; i++) {
            t_value_axis x = t_values[i];

            int16_t j = i;
            while (j > 0 && t_values[j-1].whichAxis > x.whichAxis)
            {
                t_values[j] = (t_value_axis){t_values[j-1].axisT,t_values[j-1].whichAxis,t_values[j-1].ID};
                j--;
            }
            t_values[j] = (t_value_axis){x.axisT,x.whichAxis,x.ID};
        }

        for(int16_t i = 1; i < 5; i++) {
            t_value_axis x = t_values[i];

            int16_t j = i;
            while (j > 0 && !(!(t_values[j-1].axisT > x.axisT)))
            {
                t_values[j] = (t_value_axis){t_values[j-1].axisT,t_values[j-1].whichAxis,t_values[j-1].ID};
                j--;
            }
            t_values[j] = (t_value_axis){x.axisT,x.whichAxis,x.ID};
        }
        

        int16_t startIndex = 0;

        for (int16_t i = 0; i < 5; i++)
        {
            if(t_values[i].ID == 2) {
                startIndex = i+1;
                break;
            }
            switch (t_values[i].ID)
            {
            case 0:
                inital_position[0] += t_values[i].axisT >= 0.0f ? ray.signDir[0] : 0;
                break;
            case 1:
                inital_position[1] += t_values[i].axisT >= 0.0f ? ray.signDir[1] : 0;
                break;
            default:
                break;
            }
        }

        // printf("%u \n", current_level);
        // printf("%f \n", t.axisT);
        // printf("%f \n", startPosition.endT);
        // print_t_values(t_values);

        uint8_t noOccupied = 1;

        for (int16_t i = startIndex; i < 5; i++)
        {
            //
            // printf("%f, %f\n", xplane, yplane);
            // printf("%f, %f\n", inital_position[0], inital_position[1]);
            //printf("%f, %f\n", cellPosition[0], cellPosition[1]);
            uint8_t occupied = readGrid(cellPosition[0] + inital_position[0],cellPosition[1] + inital_position[1],current_level);
            if(occupied == 1) {
                current_level--;
                
                if(current_level == -1) {
                    
                    output.hit = OCCUPIED;
                    output.t = t.axisT;
                    output.endposition[0] = cellPosition[0] + inital_position[0];
                    output.endposition[1] = cellPosition[1] + inital_position[1];
                    if(t.whichAxis == 0) {
                        output.normal[0] = -ray.signDir[0];
                        output.normal[1] = 0.0f;
                    }
                    else {
                        output.normal[0] = 0.0f;
                        output.normal[1] = -ray.signDir[1];
                    }
                    return output;
                }
                cellPosition[0] = (cellPosition[0] + inital_position[0]) * 2.0f;
                cellPosition[1] = (cellPosition[1] + inital_position[1]) * 2.0f;
                if(inital_position[0] == 1.0) {
                    current_lower_bounds[0] = xplane;
                }
                else {
                    current_upper_bounds[0] = xplane;
                }
                if(inital_position[1] == 1.0) {
                    current_lower_bounds[1] = yplane;
                }
                else {
                    current_upper_bounds[1] = yplane;
                }
                noOccupied = 0;
                break;
            }
            else if(occupied == 0) {
                t.axisT = t_values[i].axisT;
                t.whichAxis = t_values[i].whichAxis;

                switch (t_values[i].whichAxis)
                {
                case 0:
                    inital_position[0] += ray.signDir[0];
                    break;
                case 1: 
                    inital_position[1] += ray.signDir[1];
                    break;
                default:
                    break; 
                }  
            }   
            else {
                printf(" ----- OUTSIDE OUTSIDE OUTSIDE ----- \n");
            }

            if(cellPosition[0] < 0 || cellPosition[0] >= width[current_level] || cellPosition[1] < 0 || cellPosition[1] >= height[current_level]) {
                return output;
            }

            if(t_values[i].ID == 3) {
                break;
            }
        }

        if(noOccupied == 1) {
            current_level = total_amount_of_levels-1;
            cellPosition[0] = 0.0f;
            cellPosition[1] = 0.0f;
            current_lower_bounds[0] = 0.0f;
            current_lower_bounds[1] = 0.0f;
            current_upper_bounds[0] = start_upper_bounds[0];
            current_upper_bounds[1] = start_upper_bounds[1];
        }
    }

    return output;

}

TraversalOutput hierarchicalPlanes4Total(RayInput ray, AABBIntersectOutput startPosition) {
    TraversalOutput result = {0};

    float t = startPosition.startT;
    float cellPosition[2] = {startPosition.cellPosition[0],startPosition.cellPosition[1]};
    float totalT[2] = {0};
    float cellSize = 1.0f;
    int32_t currentLevel = 0;
    int32_t mask[2] = {0,0};
    float planes[2] = {0};

    if(startPosition.enteringDirection[0] == 1.0f) {
        mask[0] = 1.0f;
    }
    else {
        mask[1] = 1.0f;
    }

    while(cellPosition[0] >= 0.0f && cellPosition[0] < (float)width[0] 
    && cellPosition[1] >= 0.0f && cellPosition[1] < (float)height[0])
    {
        //printf("%f, %f\n", cellPosition[0], cellPosition[1]);

        currentLevel = readOctree(cellPosition);

        if(currentLevel == -1) {
            result.endposition[0] = cellPosition[0];
            result.endposition[1] = cellPosition[1];
            result.hit = OCCUPIED;
            result.normal[0] = -((float)mask[0]) * ray.signDir[0];
            result.normal[1] = -((float)mask[1]) * ray.signDir[1];
            result.t = t;
            break;
        }

        cellSize = ((float)(1 << currentLevel));

        planes[0] = (floorf(cellPosition[0]/cellSize) + ray.isDirPos[0]) * cellSize;
        planes[1] = (floorf(cellPosition[1]/cellSize) + ray.isDirPos[1]) * cellSize;

        totalT[0] = (planes[0] - ray.rayPosition[0]) * ray.invDir[0];
        totalT[1] = (planes[1] - ray.rayPosition[1]) * ray.invDir[1];

        mask[0] = totalT[0] <= totalT[1] ? 1 : 0;
        mask[1] = totalT[1] < totalT[0] ? 1 : 0;

        if (mask[0])
        {
            t = totalT[0];
            cellPosition[0] = planes[0] - ray.isDirNeg[0];

            cellPosition[1] = 
            ray.signDir[1] * fminf
            (
            ray.signDir[1]*planes[1]-ray.isDirPos[1],
                fmaxf
                (
                    0.0f,
                    ray.signDir[1]*
                    (
                        floorf
                        (
                            ray.rayDir[1]*t+ray.rayPosition[1]
                        )
                        -cellPosition[1]
                    )
                )
                +ray.signDir[1]*cellPosition[1]
            );
        }
        else {
            t = totalT[1];
            cellPosition[1] = planes[1] - ray.isDirNeg[1];

            cellPosition[0] = ray.signDir[0]*fminf(
            ray.signDir[0]*planes[0]-ray.isDirPos[0],
            fmaxf(
                0.0f,
                ray.signDir[0]*(
                    floorf(ray.rayDir[0]*t+ray.rayPosition[0])-cellPosition[0]
                )
            )
            +ray.signDir[0]*cellPosition[0]);
           
        }

    } 

    return result;
}

TraversalOutput hierarchicalPlanes8TotalPost(RayInput ray, AABBIntersectOutput startPosition) {
    TraversalOutput result = {0};

    float t = startPosition.startT;
    float cellPosition[2] = {startPosition.cellPosition[0],startPosition.cellPosition[1]};
    float totalT[2] = {0};
    float cellSize = 1.0f;
    int32_t currentLevel = 0;
    int32_t mask[2] = {0,0};
    float planes[2] = {0};

    if(startPosition.enteringDirection[0] == 1.0f) {
        mask[0] = 1.0f;
    }
    else {
        mask[1] = 1.0f;
    }

    while(cellPosition[0] >= 0.0f && cellPosition[0] < (float)width[0] 
    && cellPosition[1] >= 0.0f && cellPosition[1] < (float)height[0])
    {
        //printf("%f, %f\n", cellPosition[0], cellPosition[1]);

        currentLevel = readOctree(cellPosition);

        if(currentLevel == -1) {
            result.endposition[0] = floorf(cellPosition[0]);
            result.endposition[1] = floorf(cellPosition[1]);
            result.hit = OCCUPIED;
            result.normal[0] = -((float)mask[0]) * ray.signDir[0];
            result.normal[1] = -((float)mask[1]) * ray.signDir[1];
            result.t = t;
            break;
        }

        cellSize = ((float)(1 << currentLevel));

        planes[0] = (floorf(cellPosition[0]/cellSize) + ray.isDirPos[0]) * cellSize;
        planes[1] = (floorf(cellPosition[1]/cellSize) + ray.isDirPos[1]) * cellSize;

        totalT[0] = (planes[0] - ray.rayPosition[0]) * ray.invDir[0];
        totalT[1] = (planes[1] - ray.rayPosition[1]) * ray.invDir[1];

        mask[0] = totalT[0] <= totalT[1] ? 1 : 0;
        mask[1] = totalT[1] < totalT[0] ? 1 : 0;

    

        if(mask[0]) {
            cellPosition[1] += 
            ray.rayDir[1]*fmaxf(0.0,totalT[0] - (cellPosition[1] - ray.rayPosition[1]) * ray.invDir[1]);
            t = totalT[0];
            cellPosition[0] = planes[0] - ray.isDirNeg[0];
            
        }
        else {
            cellPosition[0] += 
            ray.rayDir[0]*fmaxf(0.0,totalT[1] - (cellPosition[0] - ray.rayPosition[0]) * ray.invDir[0]);
            t = totalT[1];
            cellPosition[1] = planes[1] - ray.isDirNeg[1];
        }
    } 

    return result;
}

TraversalOutput normalDDA(RayInput ray, AABBIntersectOutput startPosition)
{
    TraversalOutput result = {0};

    float t = startPosition.startT;
    float cellPosition[2] = {startPosition.cellPosition[0],startPosition.cellPosition[1]};

    float totalT[2] = {(cellPosition[0] + ray.isDirPos[0] - ray.rayPosition[0]) * ray.invDir[0], 
                       (cellPosition[1] + ray.isDirPos[1] - ray.rayPosition[1]) * ray.invDir[1]};

    int32_t mask[2] = {0};

    if(startPosition.enteringDirection[0] == 1.0f) {
        mask[0] = 1.0f;
    }
    else {
        mask[1] = 1.0f;
    }

    while(cellPosition[0] >= 0.0f && cellPosition[0] < (float)width[0] 
    && cellPosition[1] >= 0.0f && cellPosition[1] < (float)height[0]) {
        uint8_t occupiedCell = readGrid((int32_t)cellPosition[0],(int32_t)cellPosition[1], 0);

        //printf("Normal DDA pos: %.0f %.0f\n", cellPosition[0], cellPosition[1]);

        if(occupiedCell) {
            
            result.endposition[0] = cellPosition[0];
            result.endposition[1] = cellPosition[1];
            result.hit = OCCUPIED;
            result.normal[0] = -((float)mask[0]) * ray.signDir[0];
            result.normal[1] = -((float)mask[1]) * ray.signDir[1];
            result.t = t;
            break;
        }

        totalT[0] = (cellPosition[0] + ray.isDirPos[0] - ray.rayPosition[0]) * ray.invDir[0];
        totalT[1] = (cellPosition[1] + ray.isDirPos[1] - ray.rayPosition[1]) * ray.invDir[1];

        mask[0] = totalT[0] <= totalT[1] ? 1 : 0;
        mask[1] = totalT[1] < totalT[0] ? 1 : 0;

        if(totalT[0] <= totalT[1]) {
            cellPosition[0] += ray.signDir[0];
            t = totalT[0];
            //totalT[0] += ray.absDir[0];
        }
        else {
            cellPosition[1] += ray.signDir[1];
            t = totalT[1];
            //totalT[1] += ray.absDir[1];
        }

    }

    return result;
}

AABBIntersectOutput aabbtest(uint32_t startWidth, uint32_t startHeight,
                             RayInput ray)
{
    AABBIntersectOutput result = {0};

    uint8_t success = 1;
    float t0 = 0.0f;
    float t1 = INFINITY;

    float lowerLeftCorner[2] = {0.0f, 0.0f};
    float upperRightCorner[2] = {(float)startWidth, (float)startHeight};

    float tNear = (lowerLeftCorner[0] - ray.rayPosition[0]) * ray.invDir[0];
    float tFar = (upperRightCorner[0] - ray.rayPosition[0]) * ray.invDir[0];

    // tNear > tFar
    if (ray.signDir[0] == -1.0f)
    {
        float temp = tFar;
        tFar = tNear;
        tNear = temp;
    }

    if (tNear > t0)
    {
        t0 = tNear;
        result.enteringDirection[0] = 1.0f;
        result.enteringDirection[1] = 0.0f;
    }
    // t0 = tNear > t0 ? tNear : t0;

    if (tFar < t1)
    {
        t1 = tFar;
        result.exitingDirection[0] = 1.0f;
        result.exitingDirection[1] = 0.0f;
    }
    // t1 = tFar < t1 ? tFar : t1;

    success = t0 > t1 ? 0 : success;

    tNear = (lowerLeftCorner[1] - ray.rayPosition[1]) * ray.invDir[1];
    tFar = (upperRightCorner[1] - ray.rayPosition[1]) * ray.invDir[1];

    if (ray.signDir[1] == -1.0f)
    {
        float temp = tFar;
        tFar = tNear;
        tNear = temp;
    }

    if (tNear > t0)
    {
        t0 = tNear;
        result.enteringDirection[0] = 0.0f;
        result.enteringDirection[1] = 1.0f;
    }
    // t0 = tNear > t0 ? tNear : t0;
    if (tFar < t1)
    {
        t1 = tFar;
        result.exitingDirection[0] = 0.0f;
        result.exitingDirection[1] = 1.0f;
    }
    // t1 = tFar < t1 ? tFar : t1;

    success = t0 > t1 ? 0 : success;

    if (success)
    {
        if(result.enteringDirection[0] == 1.0f) {
            
            result.cellPosition[0] = ray.isDirNeg[0] == 1.0f ? upperRightCorner[0] - 1.0f : lowerLeftCorner[0]; 

            float yposition = 0.0f;

            int32_t level = total_amount_of_levels-1;

            float ylowerbound = (float)(lowerLeftCorner[1]);
            float yupperbound = (float)(upperRightCorner[1]);

            if(ray.isDirNeg[1] == 1.0) {
                ylowerbound = (float)(upperRightCorner[1]);
                yupperbound = (float)(lowerLeftCorner[1]);
            }

            while(level >= 0) {
                float m = (ylowerbound + yupperbound) / 2.0f;
                float whatToAddNext = (m-ray.rayPosition[1]) <= 0.0f ? 1.0f : 0.0f;
                float mray = (m-ray.rayPosition[1])*ray.invDir[1];
                if(m-ray.rayPosition[1] == 0.0f) {
                    whatToAddNext = ray.isDirPos[1];
                }
                if(mray < t0 && mray > 0.0f) {
                    whatToAddNext += ray.signDir[1];
                }
                if(whatToAddNext == 1.0f) {
                    if(ray.isDirNeg[1] == 1.0) {
                        yupperbound = m;
                    }
                    else {
                        ylowerbound = m;
                    }
                }
                else {
                    if(ray.isDirNeg[1] == 1.0) {
                        ylowerbound = m;
                    }
                    else {
                        yupperbound = m;
                    }
                }
                yposition = yposition * 2 + whatToAddNext;
                level--;
            }
            
            result.cellPosition[1] = yposition;
            
        }
        else {
            result.cellPosition[1] = ray.isDirNeg[1] == 1.0f ? upperRightCorner[1] - 1.0f : lowerLeftCorner[1]; 

            float xposition = 0.0f;

            int32_t level = total_amount_of_levels-1;

            float xlowerbound = (float)(lowerLeftCorner[0]);
            float xupperbound = (float)(upperRightCorner[0]);

            if(ray.isDirNeg[0] == 1.0) {
                xlowerbound = (float)(upperRightCorner[0]);
                xupperbound = (float)(lowerLeftCorner[0]);
            }

            while(level >= 0) {
                float m = (xlowerbound + xupperbound) / 2.0f;
                float whatToAddNext = (m-ray.rayPosition[0]) <= 0.0f ? 1.0f : 0.0f;
                float mray = (m-ray.rayPosition[0])*ray.invDir[0];
                if(m-ray.rayPosition[0] == 0.0f) {
                    whatToAddNext = ray.isDirPos[0];
                }
                if(mray <= t0 && mray > 0.0f) {
                    whatToAddNext += ray.signDir[0];
                }
                if(whatToAddNext == 1.0f) {
                    if(ray.isDirNeg[0]) {
                        xupperbound = m;
                    } 
                    else
                    {
                        xlowerbound = m;
                    }
                }
                else {
                    if(ray.isDirNeg[0]) {
                        xlowerbound = m;
                    } 
                    else
                    {
                        xupperbound = m;
                    }
                    
                }
                xposition = xposition * 2 + whatToAddNext;
                level--;
            }
            
            result.cellPosition[0] = xposition;

            
        }

        //printf("here: %f %f\n", result.cellPosition[0], result.cellPosition[1]);

        result.startT = t0;
        result.endT = t1;
        result.hitOrNot = 1;
    }
    else
    {
        result.hitOrNot = 0;
    }

    return result;
}

void printOctree(uint32_t level) {
    for (int32_t i = height[level]-1u; i >= 0; i--)
    {
        for (uint32_t j = 0; j < width[level]; j++)
        {
            printf("|%u", grid[level][j + width[level]*i]);
        }
        printf("|\n");
    }
    printf("\n");
}



void printScore(MethodScore score, char methodName[]) {
    printf("%s score:\n", methodName);
    printf("Position %llu/%llu\n", score.correctPositionChecks, score.amountPositionChecks);
    printf("Normal %llu/%llu\n", score.correctNormalChecks, score.amountNormalChecks);
}

void writeToGrid(uint16_t grid_number) {
    if(alreadyAllocated == 0) {
        uint32_t startWidth = 4;
        uint32_t startHeight = 4;

        total_amount_of_levels = 2;

        height = malloc(sizeof(uint32_t)*total_amount_of_levels);
        width = malloc(sizeof(uint32_t)*total_amount_of_levels);

        for (uint32_t i = 0; i < total_amount_of_levels; i++)
        {
            width[i] = startWidth;
            height[i] = startHeight;

            //printf("%u %u\n", width[i], height[i]);

            startWidth =  startWidth  - (startWidth >> 1u);
            startHeight = startHeight - (startHeight >> 1u);
        }


        grid = malloc(sizeof(uint8_t*)*total_amount_of_levels);

        for (uint32_t i = 0; i < total_amount_of_levels; i++)
        {
            grid[i] = malloc(sizeof(uint8_t)*height[i]*width[i]);
        }
    }

    for (uint64_t y = 0; y < 4; y++)
    {
        for (uint64_t x = 0; x < 4; x++)
        {
            grid[0][x + width[0] * y] = (grid_number & 1u);
            grid_number = grid_number >> 1u;
        }
        
    }

    for (uint32_t i = 1; i < total_amount_of_levels; i++)
    {
        makeNewOctreeLevel(i);
    }

    alreadyAllocated = 1;
}

void randomize_grid() {
    if(alreadyAllocated) {
        free(height);
        free(width);
        for (uint32_t i = 0; i < total_amount_of_levels; i++)
        {
            free(grid[i]);
        }
        free(grid);
    }

    uint32_t startWidth = (1 << (random() % 8 + 1));//random() % 64) + 4;
    uint32_t startHeight = startWidth; //(random() % 64) + 4;
            
    total_amount_of_levels = ((uint32_t)ceilf(log2f((float)__max(startWidth,startHeight))));

    height = malloc(sizeof(uint32_t)*total_amount_of_levels);
    width = malloc(sizeof(uint32_t)*total_amount_of_levels);

    for (uint32_t i = 0; i < total_amount_of_levels; i++)
    {
        width[i] = startWidth;
        height[i] = startHeight;

        //printf("%u %u\n", width[i], height[i]);

        startWidth =  startWidth  - (startWidth >> 1u);
        startHeight = startHeight - (startHeight >> 1u);
    }
    

    grid = malloc(sizeof(uint8_t*)*total_amount_of_levels);

    for (uint32_t i = 0; i < total_amount_of_levels; i++)
    {
        grid[i] = malloc(sizeof(uint8_t)*height[i]*width[i]);
    }

    for (uint16_t i = 0; i < height[0]; i++)
    {
        for (uint16_t j = 0; j < width[0]; j++)
        {
            grid[0][j + i * width[0]] = (random() % 8) / 7;
        }
    }

    for (uint32_t i = 1; i < total_amount_of_levels; i++)
    {
        makeNewOctreeLevel(i);
    }
    

    alreadyAllocated = 1;
}

void randomize_ray_position(RayInput* rayInput) {
    float rayPosition[2] = {0.0f, 0.0f};
    
    float zeroToOne = ((float)random()/4294967296.0f);

    /*
       /0\
       1 2
       \3/
    */

    switch(random() % 4) {
        case 0:
            rayPosition[0] = ((float)(width[0] + 2) * zeroToOne) - 1.0f;
            rayPosition[1] = (float)(height[0] + 1);      
            break;
        case 1: 
            rayPosition[0] = -1.0f;
            rayPosition[1] = ((float)(height[0] + 2) * zeroToOne) - 1.0f;
            break;
        case 2:
            rayPosition[0] = (float)(width[0] + 1);
            rayPosition[1] = ((float)(height[0] + 2) * zeroToOne) - 1.0f;
            break;
        case 3:
            rayPosition[0] = ((float)(width[0] + 2) * zeroToOne) - 1.0f;
            rayPosition[1] = -1.0f;
            break;
    }

    // x-direction
    rayInput->rayPosition[0] = rayPosition[0];
    
    // y-direction
    rayInput->rayPosition[1] = rayPosition[1];

}

void randomize_ray_direction(RayInput* rayInput) {
    
    float rayDir[2] = {0.0f, 0.0f};      // randomize and normalize

    float randomPosInGrid[2] = {0.0f,0.0f};

    float zeroToOne = ((float)random()/4294967296.0f);

    randomPosInGrid[0] = zeroToOne * (float)width[0];
    zeroToOne = ((float)random()/4294967296.0f);
    randomPosInGrid[1] = zeroToOne * (float)height[0];

    rayDir[0] = randomPosInGrid[0] - rayInput->rayPosition[0];
    rayDir[1] = randomPosInGrid[1] - rayInput->rayPosition[1];

    float length = sqrtf(((rayDir[0] * rayDir[0]) + (rayDir[1] * rayDir[1])));
    rayDir[0] /= length;
    rayDir[1] /= length;


    float invDir[2] = {1.0f / rayDir[0], 1.0f / rayDir[1]};
    float absDir[2] = {fabsf(invDir[0]), fabsf(invDir[1])};
    float isDirPos[2] = {rayDir[0] >= 0.0f ? 1.0f : 0.0f, rayDir[1] >= 0.0f ? 1.0f : 0.0f};
    float isDirNeg[2] = {rayDir[0] >= 0.0f ? 0.0f : 1.0f, rayDir[1] >= 0.0f ? 0.0f : 1.0f};
    float isDirNegAlt[2] = {rayDir[0] > 0.0f ? 0.0f : 1.0f, rayDir[1] > 0.0f ? 0.0f : 1.0f};
    float signDir[2] = {rayDir[0] >= 0.0f ? 1.0f : -1.0f, rayDir[1] >= 0.0f ? 1.0f : -1.0f};
    
    
    // x-direction
    rayInput->rayDir[0] = rayDir[0];
    rayInput->invDir[0] = invDir[0];
    rayInput->absDir[0] = absDir[0];
    rayInput->isDirPos[0] = isDirPos[0];
    rayInput->isDirNeg[0] = isDirNeg[0];
    rayInput->isDirNegAlt[0] = isDirNegAlt[0];
    rayInput->signDir[0] = signDir[0];
    
    // y-direction
    rayInput->rayDir[1] = rayDir[1];
    rayInput->invDir[1] = invDir[1];
    rayInput->absDir[1] = absDir[1];
    rayInput->isDirPos[1] = isDirPos[1];
    rayInput->isDirNeg[1] = isDirNeg[1];
    rayInput->isDirNegAlt[1] = isDirNegAlt[1];
    rayInput->signDir[1] = signDir[1];

    //printf("ray pos %.30f %.30f\n", rayInput->rayPosition[0], rayInput->rayPosition[1]);

    //printf("direction %.30f %.30f\n", rayInput->rayDir[0], rayInput->rayDir[1]);
}

void specific_ray_direction(RayInput* rayInput) {
    
    float rayDir[2] = {rayInput->rayDir[0], rayInput->rayDir[1]};

    float length = sqrtf(((rayDir[0] * rayDir[0]) + (rayDir[1] * rayDir[1])));
    rayDir[0] /= length;
    rayDir[1] /= length;


    float invDir[2] = {1.0f / rayDir[0], 1.0f / rayDir[1]};
    float absDir[2] = {fabsf(invDir[0]), fabsf(invDir[1])};
    float isDirPos[2] = {rayDir[0] >= 0.0f ? 1.0f : 0.0f, rayDir[1] >= 0.0f ? 1.0f : 0.0f};
    float isDirNeg[2] = {rayDir[0] >= 0.0f ? 0.0f : 1.0f, rayDir[1] >= 0.0f ? 0.0f : 1.0f};
    float isDirNegAlt[2] = {rayDir[0] > 0.0f ? 0.0f : 1.0f, rayDir[1] > 0.0f ? 0.0f : 1.0f};
    float signDir[2] = {rayDir[0] >= 0.0f ? 1.0f : -1.0f, rayDir[1] >= 0.0f ? 1.0f : -1.0f};
    
    
    // x-direction
    rayInput->rayDir[0] = rayDir[0];
    rayInput->invDir[0] = invDir[0];
    rayInput->absDir[0] = absDir[0];
    rayInput->isDirPos[0] = isDirPos[0];
    rayInput->isDirNeg[0] = isDirNeg[0];
    rayInput->isDirNegAlt[0] = isDirNegAlt[0];
    rayInput->signDir[0] = signDir[0];
    
    // y-direction
    rayInput->rayDir[1] = rayDir[1];
    rayInput->invDir[1] = invDir[1];
    rayInput->absDir[1] = absDir[1];
    rayInput->isDirPos[1] = isDirPos[1];
    rayInput->isDirNeg[1] = isDirNeg[1];
    rayInput->isDirNegAlt[1] = isDirNegAlt[1];
    rayInput->signDir[1] = signDir[1];

    //printf("ray pos %.30f %.30f\n", rayInput->rayPosition[0], rayInput->rayPosition[1]);

    //printf("direction %.30f %.30f\n", rayInput->rayDir[0], rayInput->rayDir[1]);
}

RayInput directionFromAngle(float angle) {
    RayInput rayInput;
    float rayDir[2] = {cosf(angle), sinf(angle)}; 

    float invDir[2] = {1.0f / rayDir[0], 1.0f / rayDir[1]};
    float absDir[2] = {fabsf(invDir[0]), fabsf(invDir[1])};
    float isDirPos[2] = {rayDir[0] >= 0.0f ? 1.0f : 0.0f, rayDir[1] >= 0.0f ? 1.0f : 0.0f};
    float isDirNeg[2] = {rayDir[0] >= 0.0f ? 0.0f : 1.0f, rayDir[1] >= 0.0f ? 0.0f : 1.0f};
    float signDir[2] = {rayDir[0] >= 0.0f ? 1.0f : -1.0f, rayDir[1] >= 0.0f ? 1.0f : -1.0f};
    float isDirNegAlt[2] = {rayDir[0] > 0.0f ? 0.0f : 1.0f, rayDir[1] > 0.0f ? 0.0f : 1.0f};
    
    // x-direction
    rayInput.rayDir[0] = rayDir[0];
    rayInput.invDir[0] = invDir[0];
    rayInput.absDir[0] = absDir[0];
    rayInput.isDirPos[0] = isDirPos[0];
    rayInput.isDirNeg[0] = isDirNeg[0];
    rayInput.isDirNegAlt[0] = isDirNegAlt[0];
    rayInput.signDir[0] = signDir[0];
    
    // y-direction
    rayInput.rayDir[1] = rayDir[1];
    rayInput.invDir[1] = invDir[1];
    rayInput.absDir[1] = absDir[1];
    rayInput.isDirPos[1] = isDirPos[1];
    rayInput.isDirNeg[1] = isDirNeg[1];
    rayInput.isDirNegAlt[1] = isDirNegAlt[1];
    rayInput.signDir[1] = signDir[1];

    return rayInput;

}

float angleFromDirection(RayInput rayInput) {
    return (float)(acos((double)rayInput.rayDir[0])*(double)rayInput.signDir[1]);
}

uint32_t rayQuadrant(RayInput ray) {
    uint32_t quadrant = ((uint32_t)ray.isDirNeg[0]) + 2u * ((uint32_t)ray.isDirNeg[1]);
    const uint32_t answer[4] = {0,1,3,2};
    return answer[quadrant];
}

uint8_t comparison_of_ray(RayInput candidateRay, RayInput mainRay, uint8_t intendedDirection)
{
    uint8_t is_left = 0;
    uint8_t is_right = 0;

    /*
    1|0
    ---
    2|3
    */

    uint32_t candQuadrant = rayQuadrant(candidateRay);
    uint32_t mainQuadrant = rayQuadrant(mainRay) * 4u;

    switch (candQuadrant + mainQuadrant)
    {
    case 0:
        //  cand 0 main 0
        //  depends
        if(candidateRay.rayDir[0] <= mainRay.rayDir[0] && candidateRay.rayDir[1] >= mainRay.rayDir[1] ) {
            is_left = 1;
        }
        if(candidateRay.rayDir[0] >= mainRay.rayDir[0] && candidateRay.rayDir[1] <= mainRay.rayDir[1] ) {
            is_right = 1;
        }
        break;
    case 1:
        //  cand 1 main 0
        is_left = 1;
        is_right = 0;
        break;
    case 2:
        //  cand 2 main 0
        //  depends
        if(-candidateRay.rayDir[0] >= mainRay.rayDir[0] && -candidateRay.rayDir[1] <= mainRay.rayDir[1] ) {
            is_left = 1;
        }
        if(-candidateRay.rayDir[0] <= mainRay.rayDir[0] && -candidateRay.rayDir[1] >= mainRay.rayDir[1] ) {
            is_right = 1;
        }
        break;
    case 3:
        //  cand 3 main 0
        is_left = 0;
        is_right = 1;
        break;
    case 4:
        //  cand 0 main 1
        is_left = 0;
        is_right = 1;
        break;
    case 5:
        //  cand 1 main 1
        //  depends
        if(candidateRay.rayDir[0] <= mainRay.rayDir[0] && candidateRay.rayDir[1] <= mainRay.rayDir[1] ) {
            is_left = 1;
        }
        if(candidateRay.rayDir[0] >= mainRay.rayDir[0] && candidateRay.rayDir[1] >= mainRay.rayDir[1] ) {
            is_right = 1;
        }
        break;
    case 6:
        //  cand 2 main 1
        is_left = 1;
        is_right = 0;
        break;
    case 7:
        //  cand 3 main 1
        //  depends
        if(-candidateRay.rayDir[0] >= mainRay.rayDir[0] && -candidateRay.rayDir[1] >= mainRay.rayDir[1] ) {
            is_left = 1;
        }
        if(-candidateRay.rayDir[0] <= mainRay.rayDir[0] && -candidateRay.rayDir[1] <= mainRay.rayDir[1] ) {
            is_right = 1;
        }
        break;
    case 8:
        //  cand 0 main 2
        //  depends
        if(-candidateRay.rayDir[0] <= mainRay.rayDir[0] && -candidateRay.rayDir[1] >= mainRay.rayDir[1] ) {
            is_left = 1;
        }
        if(-candidateRay.rayDir[0] >= mainRay.rayDir[0] && -candidateRay.rayDir[1] <= mainRay.rayDir[1] ) {
            is_right = 1;
        }
        break;
    case 9:
        //  cand 1 main 2
        is_left = 0;
        is_right = 1;
        break;
    case 10:
        //  cand 2 main 2
        //  depends
        if(candidateRay.rayDir[0] >= mainRay.rayDir[0] && candidateRay.rayDir[1] <= mainRay.rayDir[1] ) {
            is_left = 1;
        }
        if(candidateRay.rayDir[0] <= mainRay.rayDir[0] && candidateRay.rayDir[1] >= mainRay.rayDir[1] ) {
            is_right = 1;
        }
        break;
    case 11:
        //  cand 3 main 2
        is_left = 1;
        is_right = 0;
        break;
    case 12:
        //  cand 0 main 3
        is_left = 1;
        is_right = 0;
        break;
    case 13:
        //  cand 1 main 3
        //  depends
        if(-candidateRay.rayDir[0] <= mainRay.rayDir[0] && -candidateRay.rayDir[1] <= mainRay.rayDir[1] ) {
            is_left = 1;
        }
        if(-candidateRay.rayDir[0] >= mainRay.rayDir[0] && -candidateRay.rayDir[1] >= mainRay.rayDir[1] ) {
            is_right = 1;
        }
        break;
    case 14:
        //  cand 2 main 3
        is_left = 0;
        is_right = 1;
        break;
    case 15:
        //  cand 3 main 3
        //  depends
        if(candidateRay.rayDir[0] >= mainRay.rayDir[0] && candidateRay.rayDir[1] >= mainRay.rayDir[1] ) {
            is_left = 1;
        }
        if(candidateRay.rayDir[0] <= mainRay.rayDir[0] && candidateRay.rayDir[1] <= mainRay.rayDir[1] ) {
            is_right = 1;
        }
        break;
    }

    if(intendedDirection == 0) {
        is_left = is_left == is_right ? 0 : is_left;
        return is_left;
    }
    else {
        is_right = is_left == is_right ? 0 : is_right;
        return is_right;
    }
}

int main()
{
    setlocale(LC_NUMERIC, "French_Canada.1252");

    struct timespec t1, t2;
    struct timespec result1 = {0}, result2 = {0}, result3 = {0}, result4 = {0};
    struct timespec minresult1 = {0}, maxresult1 = {0};
    struct timespec minresult2 = {0}, maxresult2 = {0};
    struct timespec minresult3 = {0}, maxresult3 = {0};
    struct timespec minresult4 = {0}, maxresult4 = {0};
    struct timespec tempdiff;

    int seedsuccess = 0;
    unsigned long long randvalue = 0;

    seedsuccess = _rdrand64_step(&randvalue);

    if (!seedsuccess)
    {
        randvalue = t1.tv_nsec + t1.tv_sec;
    }

    zxorshiftstate = randvalue;

    zxorshiftstate = 16516837009037423286llu;

    printf("--- exact cases ---\n\n");

    RayInput rayInput;

    ray_info exactCases[] = {
        
        {-1.0f,0.0f,1.0f,0.0f},
        {-1.0f,1.0f,1.0f,0.0f},
        {-1.0f,2.0f,1.0f,0.0f},
        {-1.0f,3.0f,1.0f,0.0f},

        {-1.0f,-1.0f,1.0f,1.0f},

        {-1.0f,0.0f,1.0f,1.0f},
        {-1.0f,1.0f,1.0f,1.0f},
        {-1.0f,2.0f,1.0f,1.0f},

        {0.0f,-1.0f,0.0f,1.0f},
        {1.0f,-1.0f,0.0f,1.0f},
        {2.0f,-1.0f,0.0f,1.0f},
        {3.0f,-1.0f,0.0f,1.0f},

        {0.0f,-1.0f,1.0f,1.0f},
        {1.0f,-1.0f,1.0f,1.0f},
        {2.0f,-1.0f,1.0f,1.0f},

        {5.0f,0.0f,-1.0f,0.0f},
        {5.0f,1.0f,-1.0f,0.0f},
        {5.0f,2.0f,-1.0f,0.0f},
        {5.0f,3.0f,-1.0f,0.0f},

        {5.0f,-1.0f,-1.0f,1.0f},

        {5.0f,0.0f,-1.0f,1.0f},
        {5.0f,1.0f,-1.0f,1.0f},
        {5.0f,2.0f,-1.0f,1.0f},

        {0.0f,5.0f,0.0f,-1.0f},
        {1.0f,5.0f,0.0f,-1.0f},
        {2.0f,5.0f,0.0f,-1.0f},
        {3.0f,5.0f,0.0f,-1.0f},

        {-1.0f,5.0f,1.0f,-1.0f},
        {0.0f,5.0f,1.0f,-1.0f},
        {1.0f,5.0f,1.0f,-1.0f},
        {2.0f,5.0f,1.0f,-1.0f},

        {5.0f,5.0f,-1.0f,-1.0f},

        {5.0f,4.0f,-1.0f,-1.0f},
        {5.0f,3.0f,-1.0f,-1.0f},
        {5.0f,2.0f,-1.0f,-1.0f},

        {4.0f,5.0f,-1.0f,-1.0f},
        {3.0f,5.0f,-1.0f,-1.0f},
        {2.0f,5.0f,-1.0f,-1.0f},

        {2.0f,-1.0f,-1.0f,1.0f},
        {3.0f,-1.0f,-1.0f,1.0f},
        {4.0f,-1.0f,-1.0f,1.0f},

        {-1.0f,4.0f,1.0f,-1.0f},
        {-1.0f,3.0f,1.0f,-1.0f},
        {-1.0f,2.0f,1.0f,-1.0f},
        
    };

    MethodScore HP4T_score_specific_cases = {0};
    MethodScore HP8T_score_specific_cases = {0};
    MethodScore H4T_score_specific_cases = {0};

    for (uint32_t i = 1; i < 65535; i++)
    {
        writeToGrid(i);

        for (uint32_t j = 0; j < 44; j++)
        {   
            rayInput.rayDir[0] = exactCases[j].dir_x;
            rayInput.rayDir[1] = exactCases[j].dir_y;
            rayInput.rayPosition[0] = exactCases[j].pos_x;
            rayInput.rayPosition[1] = exactCases[j].pos_y;

            specific_ray_direction(&rayInput);

            AABBIntersectOutput startPosition = aabbtest(width[0], height[0], rayInput);

            //printf("%d\n", startPosition.hitOrNot);

            if (startPosition.hitOrNot)
            {

                TraversalOutput endPointDDA = normalDDA(rayInput, startPosition);

                TraversalOutput endPointH4T = hero4Total(rayInput,startPosition);
                uint8_t H4Tfailure = compareTraversalOutput(endPointDDA, endPointH4T, &H4T_score_specific_cases);

                TraversalOutput endPointHP4T = hierarchicalPlanes4Total(rayInput, startPosition);
                uint8_t HP4Tfailure = compareTraversalOutput(endPointDDA, endPointHP4T, &HP4T_score_specific_cases);

                TraversalOutput endPointHP8T = hierarchicalPlanes8TotalPost(rayInput, startPosition);
                uint8_t HP8Tfailure = compareTraversalOutput(endPointDDA, endPointHP8T, &HP8T_score_specific_cases);

                
            }
        }
        
    }
    
    printScore(HP4T_score_specific_cases, "HP4T specific cases");
    printScore(HP8T_score_specific_cases, "HP8T specific cases");
    printScore(H4T_score_specific_cases, "H4T specific cases");

    printf("--- random cases ---\n\n");

    printf("Seed: %llu\n\n", zxorshiftstate);

    uint64_t howManyDontHitInsideGrid = 0;

    MethodScore HP4T_score = {0};
    MethodScore HP8T_score = {0};
    MethodScore H4T_score = {0};

    for (uint32_t i = 0; i < 100; i++)
    {
        randomize_grid();

        //printf("---------- I: %03u ----------\n", i);

        for (uint32_t j = 0; j < 100; j++)
        {
            randomize_ray_position(&rayInput);
            randomize_ray_direction(&rayInput);

            AABBIntersectOutput startPosition = aabbtest(width[0], height[0], rayInput);

            //printf("%d\n", startPosition.hitOrNot);

            if (startPosition.hitOrNot)
            {
                //printf("---------- J: %03u ----------\n", j);

                clock_gettime(CLOCK_MONOTONIC, &t1);

                TraversalOutput endPointDDA = normalDDA(rayInput, startPosition);

                clock_gettime(CLOCK_MONOTONIC, &t2);

                if(endPointDDA.hit == 0.0f) {
                    howManyDontHitInsideGrid++;
                }

                timespec_subtract(&tempdiff,&t2,&t1);

                timespec_min_special(&tempdiff,&minresult4, i, j);
                timespec_max_special(&tempdiff,&maxresult4, i, j);

                timespec_add(&result4,&tempdiff); 

                clock_gettime(CLOCK_MONOTONIC, &t1);

                TraversalOutput endPointHP4T = hierarchicalPlanes4Total(rayInput, startPosition);

                clock_gettime(CLOCK_MONOTONIC, &t2);

                timespec_subtract(&tempdiff,&t2,&t1);

                //printf("%lf\n", (double)tempdiff.tv_nsec / 1000.0);

                timespec_min_special(&tempdiff,&minresult1,i,j);
                timespec_max_special(&tempdiff,&maxresult1,i,j);

                timespec_add(&result1,&tempdiff); 

                uint8_t HP4Tfailure = compareTraversalOutput(endPointDDA, endPointHP4T, &HP4T_score);

                clock_gettime(CLOCK_MONOTONIC, &t1);

                TraversalOutput endPointHP8T = hierarchicalPlanes8TotalPost(rayInput, startPosition);

                clock_gettime(CLOCK_MONOTONIC, &t2);

                timespec_subtract(&tempdiff,&t2,&t1);

                //printf("%lf\n", (double)tempdiff.tv_nsec / 1000.0);

                timespec_min_special(&tempdiff,&minresult2,i,j);
                timespec_max_special(&tempdiff,&maxresult2,i,j);

                timespec_add(&result2,&tempdiff); 

                uint8_t H8TPfailure = compareTraversalOutput(endPointDDA, endPointHP8T, &HP8T_score);

                clock_gettime(CLOCK_MONOTONIC, &t1);

                TraversalOutput endPointH4T = hero4Total(rayInput,startPosition);

                clock_gettime(CLOCK_MONOTONIC, &t2);

                timespec_subtract(&tempdiff,&t2,&t1);

                //printf("%lf\n", (double)tempdiff.tv_nsec / 1000.0);

                timespec_min_special(&tempdiff,&minresult3,i,j);
                timespec_max_special(&tempdiff,&maxresult3,i,j);

                timespec_add(&result3,&tempdiff); 

                uint8_t H4Tfailure = compareTraversalOutput(endPointDDA, endPointH4T, &H4T_score);

                if(H4Tfailure == 1) {
                    printOutput(endPointH4T, "H4TH");
                    printOutput(endPointDDA, "normalDDA");
                    printf("start in grid (%.30f,%.30f)\n", startPosition.cellPosition[0], startPosition.cellPosition[1]);
                    printf("ray pos (%.30f,%.30f)\n", rayInput.rayPosition[0], rayInput.rayPosition[1]);
                    printf("ray dir (%.30f,%.30f)\n", rayInput.rayDir[0], rayInput.rayDir[1]);
                    printf("%u %u\n", width[0], height[0]);
                    printOctree(0);
                    printf("I: %u J: %u\n", i,j);
                    printf("\n\n");

                    return 0;
                }

            }
        }
    }

    printf("\n\n");
    printf("These didn't hit: %llu\n", howManyDontHitInsideGrid);


    printf("total time 1: %lld,%09ldsec taken\n", result4.tv_sec, result4.tv_nsec);
    printf("minresult1: %lld,%09ldsec taken\n", minresult4.tv_sec, minresult4.tv_nsec);
    printf("maxresult1: %lld,%09ldsec taken\n\n", maxresult4.tv_sec, maxresult4.tv_nsec);

    
    printScore(HP4T_score, "HP4T");
    printf("total time 1: %lld,%09ldsec taken\n", result1.tv_sec, result1.tv_nsec);
    printf("minresult1: %lld,%09ldsec taken\n", minresult1.tv_sec, minresult1.tv_nsec);
    printf("maxresult1: %lld,%09ldsec taken\n", maxresult1.tv_sec, maxresult1.tv_nsec);

    
    printScore(HP8T_score, "HP8T");
    printf("total time 2: %lld,%09ldsec taken\n", result2.tv_sec, result2.tv_nsec);
    printf("minresult2: %lld,%09ldsec taken\n", minresult2.tv_sec, minresult2.tv_nsec);
    printf("maxresult2: %lld,%09ldsec taken\n", maxresult2.tv_sec, maxresult2.tv_nsec);

    
    printScore(H4T_score, "H4T");
    printf("total time 3: %lld,%09ldsec taken\n", result3.tv_sec, result3.tv_nsec);
    printf("minresult3: %lld,%09ldsec taken\n", minresult3.tv_sec, minresult3.tv_nsec);
    printf("maxresult3: %lld,%09ldsec taken\n", maxresult3.tv_sec, maxresult3.tv_nsec);

    printf("\n\n");

    return 0;
}