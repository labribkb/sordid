
# Compile

This source code has been tested with C++17 and depends on the CGAL library. Experiments were run using version 5.4.1 of CGAL. On ubuntu 22.04, is available via package libcgal-dev.

Compilation requires version=>11 of g++.

To compile on Windows or Linux, type ```make```. For more detail, refer to the MakeFile.


# Run

To run the code, you must provide all the parameters expected by the code. The description of all parameters is given below. But first, here follows an example command with some default values for all parameters. 

```$ ./sordid <INPUT_LAYOUT_PATH> <ALPHA> <K> <MINIMUM_MOVEMENT> <MAX_ITER> <MAX_PASSES> <SCALE_STEP> <PRIME> <MONITOR>```

| Parameter  |  Type |  Description |
|:---:|:---:|:---:|
| ``INPUT_LAYOUT_PATH``  | string  |  path to the .txt file encoding the initial graph layout in the expected format (described below). The algorithm will output the result layout in ``INPUT_LAYOUT_PATH.forbid`` or ``INPUT_LAYOUT_PATH.forbidp`` depending on the ``PRIME`` argument |
|  ``ALPHA`` | float |  weight factor for ideal distance for both overlapped and non-overlapped pairs of nodes |
|  ``K`` | float  | additional weight factor for overlapped pairs of nodes  |
|  ``MINIMUM_MOVEMENT`` |  float |  threshold value for the optimization algorithm. In a pass in the optimization algorithm, if the sum of nodes movement is below that threshold, the pass is ended. |
|  ``MAX_ITER`` |  int |  maximum number of iterations in each pass in the optimization algorithm |
|  ``MAX_PASSES`` |  int |  maximum number of passes before existing FORBID (should not be reached)  |
|  ``SCALE_STEP`` | float  | minimal step size that stops the binary search for the optimal scale  |
|  ``PRIME`` | bool as int (0,1)  | wether to use the FORBID or FORBID' variant  |
|  ``MONITOR`` | bool as int (0,1)  | adds steps within the algorithm to record several information (such as number of overallps remaining, current stress) and save them in a text file|

Example: 

```$ ./sordid ../demo/test/dpd.txt 1 64 0.000001 30 100 0.1 0 0 ```

# Input graph format

The input graph should be a txt file in which:
* the first line contains the graph number of nodes (i.e., number of lines to parse)
* one line per node in the format ``x y <number of polygon vertices> <v_0_x> <v_0_y> ... <v_k_x> <v_k_y>`` where ``(v_i_x,v_i_y)`` is the position of the ith vertex of the node shape relative to the node position ``(x,y)``. Vertices must be in counter-clockwise order.

Example for a small graph with 4 nodes:
```
4
1 2 4 -1 -1 1 -1 1 1 -1 1
10 10 4 -1 -1 0 -1 1 1 0 -1
7 7 4 -2 -2 2 -2 2 2 -2 2
0 9 4 -10 -1 0 -1 10 1 0 -1
```

# Output graph format

The algorithm will output a file in almost the same format as its input, except that the first line is the execution time of the algorithm. The algorithm result is stored in a new file created at the same location as the input file, with the same name on which is appended the extension ``.sordid``