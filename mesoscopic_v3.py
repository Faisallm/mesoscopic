import numpy as np
import random
import math


class Mesoscale:
    
    
    def __init__(self, L, W, H, dmin, dmax, number):
        
        
        # initialization of the global background grid
        # the global background grid of concrete specimen.
        # let's take unit in mm.
        # L (length)
        # W (width)
        # H (height) of the concrete specimen.
        self.L = L
        self.W = W
        self.H = H
        
        # minimum and maximum aggregate size value.
        # this value helps in the generation of the polyhedral aggregate.
        self.dmin = dmin
        self.dmax = dmax
        
        # the mesh element size (e)
        # this is taken as (1/4) ~ (1/8) of minimum aggregate particle size dmin.
        # if the size of mesh element is too small, it will greatly increase...
        # the amount of unneccessary storage and calculation.

        # (we used the largest size possible for our dmin to decrease computation)
        # we could also use (1/8), leading to larger computation
        self.e = int((1/4)*self.dmin)
        
        # represents the number of elements along x, y and z direction.
        self.l = int(self.L//self.e)
        self.m = int(self.W//self.e)
        self.n = int(self.H//self.e)
        
        # self.number represent the number of test aggregates to be generated for placement into the new area
        self.number = number
        
        assert((self.L*self.W*self.H) == (self.l*self.e * self.m * self.e * self.n * self.e))
        
        # the total number of elements in the concrete specimen
        self.total_elements = self.l * self.m * self.n
        
        # initialization of the global background grid
        self.glob = self._global_matrix(self.l, self.m, self.n)
        
        # aggregate size ranges
        self.aggregate_size_ranges = [[5, 10], [10, 15], [15, 20]]
        
        # this will compute the cumulative probability that an aggregate..
        # passes a sieve opening with a size of d mm.
#         a = [i[1] for i in self.aggregate_size_ranges]
#         print(a)
        self.probabilities = [self._fuller_curve(i) for i in self.aggregate_size_ranges]
#         print(f"probabilities: {self.probabilities}")
        
        # generating a certain number of aggregates.
        # generate a set of aggregates with different dmin and dmax.
        # nested list comprehension.
        
        self.aggs = [self._generate_polyhedron(dmin, dmax) \
                     for dmin, dmax in random.choices([i for i in self.aggregate_size_ranges], weights=self.probabilities) \
                     for j in range(self.number)]
        
        # arranging the aggregates based on area.
        # this is done in order to achieve a higher aggregate content.
#         print(self.aggs[0])
        self.arranged_aggs = sorted(self.aggs, key=lambda x: self._polyhedral_area(x), reverse=True)
        
        # generating a certain number of aggregates.
#         self.aggs = [self._generate_polyhedron(self.dmin, self.dmax) for i in range(self.number)]
        
        self._condition = True
        self._count = 0
        
        # testing
        for i, a in enumerate(self.arranged_aggs):
            while self._condition and self._count < 8000:
                self.agg = self._translate(a,
                                          self._random_translation_point(self.l, self.m, self.n, self.e))
                # local grid
                # If their is an intrusion, I want the aggregate to be randomly translated to a new
                # location for placing for atleast 500 times.
                self.loc, self._condition = self._identify_aggregate_and_itz(self.agg)
                #print(f"Local: {self.loc}")
                # checking if the aggregate was placed
                if not self._condition:
                    # break out of loop
                    self._condition = False
                    self._count = 0
                    print("no intrusion")
                else:
                    self._count += 1;
                    print("intrusion detected!")
            # to carry out the same check for other aggregates        
            self._condition = True
            print(f"Placed #{i+1} aggregate.")
            
            
        print(f"Local: {self.loc}")
        print(f"Global: {self.glob}")
        
        
        
    def _element_coordinates(self, i: int, j: int, k: int, e):
        """This code determines the eight nodes of any element
       (i, j, k) in the global 3D background grid.
       Where
       i = element no along x-axis
       j = element no along y-axis
       k = element no along z-axis.
       The type of array returned is a numpy array to optimize
       numerical calculations.
       """

        coordinates = np.array((
            (i * e, j * e, k * e),
            ((i+1)*e, j * e, k * e),
            (i*e, (j+1)*e, k * e),
            ((i+1)*e, (j+1)*e, k * e),
            (i*e, j*e, (k+1)*e),
            ((i+1)*e, j * e, (k+1)*e),
            (i*e, (j+1)*e, (k+1)*e),
            ((i+1)*e, (j+1)*e, (k+1)*e)
        ))
        
        return coordinates
    
    def _polyhedral_area(self, vertices):
    
        if not np.allclose(vertices[0], vertices[-1]):
            vertices = np.vstack((vertices, vertices[0]))

        # Initialize area accumulator
        area = 0.0

        # Calculate the area using the Shoelace formula
        for i in range(len(vertices) - 1):
            cross_product = np.cross(vertices[i], vertices[i + 1])
            area += np.linalg.norm(cross_product)

        # Divide by 2 to get the actual area
        area /= 2.0

        return area
    
    def _fuller_curve(self, d):
        # these values do not change.
        d=d[0]
        n = 0.5
        D = 20
        return ((d/D)**n)
    
    
    def _element_central_point(self, i: int, j: int, k: int, e):
        """The coordinates of the central point of any element
        (i, j, k) in the global background grid.
        The type of array returned is a numpy array to optimize
           numerical calculations."""
        
        return np.array(((i*e) + e/2, (j*e)+e/2, (k*e)+e/2))
        
        
    def _locate_element(self, x: int, y: int, z: int, e):
        """If the coordinates of any point in the background
        mesh block are (x, y, z), the number (i, j, k) of the element
        where the point can be located is defined in this function.
        The type of array returned is a numpy array to optimize
           numerical calculations."""
        
        return np.array((int(x // e), int(y // e), int(z // e)))
    
    
    def _global_matrix(self, l: int, m: int, n: int):
        """This a 3D matrix Mijk(0 <= i <= l-1, 0 <= j <= m-1, 0<= k<= n-1).
        is used to store the material attributes of all the global
        background elements.
        It will store elements such as Mijk=1,2,3.
        Meaning that the grid element (i, j, k) is a mortar element,
        an ITZ element or an aggregate element, respectively.

        1 = mortar
        2 = ITZ (Interfacial Transition zone)
        3 = aggregate(spherical, ellipsoidal, polyhedral)

        It is an array of matrices.
        This will probably use alot of compute.
        But the idea is to produce a 3D array that you can easily
        access and change its contents.

        Initialliy all the contents of this array will be initialized 
        as mortar having a value of 1.
        """
        # initialized blank array     
        matrix = []

        for i in range(l):
            i_array = []
            for j in range(m):
                m_array = []
                for k in range(n):
                    m_array.append(1)
                i_array.append(m_array)
            matrix.append(i_array)

        return np.array(matrix)
    
    
    def _generate_points(self, radius):
        """This function is responsible for the generation of the coordinates
        (x, y, z) of the vertex of a polygon.
        radius: radius of the aggregate."""

        i_neeta = random.random()
        i_eeta = random.random()

        # azimuth angle and zenith angle
        azimuth = i_neeta * 2 * math.pi
        zenith = i_eeta * 2 * math.pi

        # coordinates
        xi = radius * math.sin(zenith) * math.cos(azimuth)
        yi = radius * math.sin(zenith) * math.sin(azimuth)
        zi = radius * math.cos(azimuth)

        return np.array((xi, yi, zi))
    
    
    def _generate_polyhedron(self, dmin, dmax):
        """
        This function will return a  single polyhedron.
        where n is the number of vertices.
        we are going to store the vertices of the coordinates in numpy.
        n: represent the number of vertices of the polygon you want to generate.
        dmin: represent the minimum aggregate size.
        dmax: represent the maximum aggregate size.

        The number of three-dimensional random points is preferably
        controlled between 15 and 25.
        """
        neeta = random.random()
        eeta = random.random()
        # radius
        r = (dmin/2) + (neeta * ((dmax-dmin)/2))
        polyhedron = []
        # as shown in the paper, the usual number of
        # polyhedral vertex fall between 15 and 25
        choices = [i for i in range(15, 25)]

        for i in range(random.choice(choices)):
            p = self._generate_points(r)
            polyhedron.append(p)

        return np.array(polyhedron)
    
    
    def _random_translation_point(self, l: int, m: int, n: int, e):
        """Represent the coordinates (xp, yp, zp) of the
        random translation points given by: (n7le, n8me, n9ne).
        where n7, n8 and n9 are random translation points
        between 0 and 1.
        l = L//e
        m = W//e
        n = H//e
        L = length of the concrete specimen
        W = Width of the concrete specimen
        H = Height of the concrete specimen
        e = element size (1/4)~(1/8) of dmin.
        To further improve the efficiency of aggregate random placement,
        the random translation point P can be controlled to fall inside the
        mortar elements.
        To ensure that the random translation point falls within the
        dimensions of the concrete specimen
        0 <= Xi"""


        while True:
            n7 = random.random()
            n8 = random.random()
            n9 = random.random()

            xp, yp, zp = (n7 * l * e, n8 * m * e, n9 * n * e)

            # ensuring that the random translation point P is within the concrete specimen.
            # if it is outside, generate another set of random translation point P.
            if not ((0 <= xp < l*e) and (0 <= yp < n * e) and (0 <= zp < m * e)):
                continue

            # get the element where this point is located
            element = self._locate_element(xp, yp, zp, e)


            # check that the element that the random translation point...
            # is located is not currently a mortar.
            # 1 indicates mortar.

            if (self.glob[element[0], element[1], element[2]] == 1):
                break

        return np.array((xp, yp, zp))
    
    
    def _translate(self, aggregate, translation_point):
        """The newly generated random polyhedral aggregate is directly put into
        the background grid space by RANDOM TRANSLATION.
        Assuming the coordinates of the random translation point are
        (xp, yp, zp)=(n7le, n8me, n9ne).
        Therefore, the coordinates of the vertices after translation are:
        Xi = xi + xp, Yi = yi + yp, Zi = zi + zp."""

        xp, yp, zp = self._random_translation_point(self.l, self.m, self.n, self.e)

        # actual process of translating the aggregate.
        # adding the translated point to each vertex of the aggregate.
        translated_aggregate = [(i[0]+xp, i[1] + yp, i[2] + zp) for i in aggregate]

        return np.array(translated_aggregate)
    
    
    def _minimum_and_maximum(self, array):
        """According to the new coordinates of the vertices of the...
        polyhedral aggregate. The minimum values (Xmin, Ymin, Zmin) and
        the maximum values (Xmax, Ymax, Zmax) of all the
        vertex coordinates in the three-dimensional direction
        are obtained in this function."""

        # axis 0 refers to the column.
        ir, jr, kr = np.amin(array, axis=0)
        Ir, Jr, Kr = np.amax(array, axis=0)

        return (int(ir), int(jr), int(kr)), (int(Ir), int(Jr), int(Kr))
    
    def _initialize_local_matrix(self, lower_limit_coordinate, upper_limit_coordinate):
        """Represents a 3D matrix B used to temporarily store
        the material attributes of the current local background grid.
        The matrix B is to be initialized to 1, that is all elements in the
        local background grid are initialized as MORTAR ELEMENTS."""

        ((x_min, y_min, z_min), (x_max, y_max, z_max)) = lower_limit_coordinate, upper_limit_coordinate

        ooo_array = []
        for i in range(x_min-1, x_max):
            oo_array = []
            for j in range(y_min-1, y_max):
                o_array = []
                for k in range(z_min-1, z_max):
                    o_array.append(1)
                oo_array.append(o_array)
            ooo_array.append(oo_array)

        return np.array(ooo_array)
    
    
    def _section_global_background_grid(self, aggregate):
        """The bounding box ((ir, jr, kr), (Ir, Jr, Kr)) of the newly placed aggregate can be
        determined using(red dashes) without considering the ITZ.

        While the coordinates ((ib, jb, kb), (Ib, Jb, Kb)) represents the 
        bounding box without considering the ITZ."""


        ((x_min, y_min, z_min), (x_max, y_max, z_max)) = self._minimum_and_maximum(aggregate)

        # element number of the minimum and maximum coordinates     
        # bounding box of the aggregate without ITZ
#         (ir, jr, kr) = (x_min//self.e, y_min//self.e, z_min//self.e)
#         (Ir, Jr, Kr) = (x_max//self.e, y_max//self.e, z_max//self.e)

        # bounding box of the newly placed aggregate (consideri)
        # I do not need this
        (ir, jr, kr) = np.max([0, x_min//self.e]), np.max([0, y_min//self.e]), np.max([0, z_min//self.e])
        (Ir, Jr, Kr) = np.min([self.l, x_max//self.e]), np.min([self.m, y_max//self.e]), np.min([self.n, z_max//self.e])

        # bounding box of the newly placed aggregate considering ITZ elements
        (ib, jb, kb) = np.max([0, ir-1]), np.max([0, jr-1]), np.max([0, kr-1])
        (Ib, Jb, Kb) = np.min([self.l, Ir+1]), np.min([self.m, Jr+1]), np.min([self.n, Kr+1])

        return ((ir, jr, kr), (Ir, Jr, Kr)), ((ib, jb, kb), (Ib, Jb, Kb)) 
    
    
    def _point_in_polyhedron(self, point, vertices):
        """
        Check if a point is inside a polyhedron defined by its vertices.

        Args:
        - point (numpy array): The point coordinates [x, y, z].
        - vertices (numpy array): The vertices of the polyhedron, where each row represents a vertex [x, y, z].

        Returns:
        - bool: True if the point is inside the polyhedron, False otherwise.
        """

        # Initialize counters for intersections with edges and faces
        intersection_count = 0

        # Iterate through each face of the polyhedron
        for i in range(len(vertices)):
            # Define the vertices of the current face
            v1 = vertices[i]
            v2 = vertices[(i + 1) % len(vertices)]

            # Check if the point is on the same side of the plane defined by the face
            # as the line segment connecting the two vertices of the face
            if (v1[1] > point[1]) != (v2[1] > point[1]):
                # Calculate the intersection point with the plane
                x_intersect = (v2[0] - v1[0]) * (point[1] - v1[1]) / (v2[1] - v1[1]) + v1[0]

                # If the intersection point is to the right of the test point, increment the intersection count
                if point[0] < x_intersect:
                    intersection_count += 1

        # If the number of intersections is odd, the point is inside the polyhedron
        return intersection_count % 2 == 1
    
    
    def _local_grid(self, lower_limit, upper_limit):
        """
        We have to generate the local background grid from the global
        background grid.

        Always remember that range excludes the outer limit, so don't
        forget to add + 1.
        """

        (ir, jr, kr), (Ir, Jr, Kr) = lower_limit, upper_limit

        ooo_array = []
        for i in range(0, (Ir-ir)+1):
            oo_array = []
            for j in range(0, (Jr-jr)+1):
                o_array = []
                for k in range(0, (Kr-kr)+1):
                    o_array.append(1)
                oo_array.append(o_array)
            ooo_array.append(oo_array)

        return np.array(ooo_array)
    
    
    def _is_adjacent_to_aggregate(self, space, i, j, k):
        """
        The purpose of this function is to identify elements
        that are adjacent to aggregate elements. Thus help
        in classifying such elements as ITZ (interfacial transition zone).
        """
    
        adjacent_positions = [
            (i+1, j, k), (i-1, j, k),   # Check neighbors along x-axis
            (i, j+1, k), (i, j-1, k),   # Check neighbors along y-axis
            (i, j, k+1), (i, j, k-1)    # Check neighbors along z-axis
        ]

        for x, y, z in adjacent_positions:
            if 0 <= x < len(space) and 0 <= y < len(space[0]) and 0 <= z < len(space[0][0]):
                if space[x][y][z] == 3:
                    return True
        return False
    
    
    def _aggregate_intrusion(self, ir: int, Ir: int, jr: int, Jr: int, kr: int, Kr: int):
        """For a newly determined aggregate element (i, j, k) in the local
        background, in order to ensure to that the new aggregate element (i, j, k)
        does not overlap or contact with the old aggregate elements, its
        corresponding element (i + ib, j + jb, k + kb) in the global
        background grid cannot be an aggregate element or ITZ element.
        
        This will return a boolean indicating whether the newly placed aggregate
        in the local background grid is interfaring with a placed old aggregate
        in the global background grid."""
        
        contains = False
        
        for i in range(ir-1, Ir):
            for j in range(jr-1, Jr):
                for k in range(kr-1, Kr):
                    if self.glob[i][j][k] == 3 or self.glob[i][j][k] == 2: 
                        contains = True
                        
        return contains
        
    
    
    def _identify_aggregate_and_itz(self, agg):
        """
        1, loop through all the elements in the local background grid of current aggregate.
        i.e all the elements in the red dotted block.

        2, According to the element number of the local background grid,
        the corresponding element of the global background grid is obtained.
        and then the coordinates of the center points are obtained (eq 4).
        """

        # a heck of unpacking took place here.
        # I had a little misunderstanding but this should be known...
        # as the subsection of an element in the globaal background grid
        ((ir, jr, kr), (Ir, Jr, Kr)), ((ib, jb, kb), (Ib, Jb, Kb)) = self._section_global_background_grid(agg)
#         print(f"ir: {ir}, Ir: {Ir}, jr: {jr}, Jr:{Jr}, kr: {kr}, Kr: {Kr}")
#         print(f"ib: {ib}, Ib: {Ib}, jb: {jb}, Jb:{Jb}, kb: {kb}, Kb: {Kb}")
        b_matrix = self._initialize_local_matrix((ir, jr, kr), (Ir, Jr, Kr))
        # real local grid
        # we generate the value of the local background grid
        # based on the global background grid.
        local_grids = self._local_grid((ib, jb, kb), (Ib, Jb, Kb))
        #print(f"local grid shape: {local_grids.shape}")
        # our aggregate intrusion test
        # the local background grid method makes...
        # checking for intrusion so simple
        intrusion = self._aggregate_intrusion(ir, Ir, jr, Jr, kr, Kr)
        intru = False

        # looping through the elements of the local background grid     
        # all elements in the red dotted block.
        for i in range(ir-1, Ir):
            for j in range(jr-1, Jr):
                for k in range(kr-1, Kr):
                    # element number of the global background grid
                    global_element = i, j, k
                    # element central point (global bakground grid)
                    center_point = self._element_central_point(i, j, k, self.e)
                    # coordinates of the local background grid
                    # relationship between the local and global grids.
                    li, lj, lk = (i-(ir-1), j-(jr-1), k-(kr-1))
                    #print(li, lj, lk)

                    # if the element is aggregate
                    if self._point_in_polyhedron(center_point, agg):
                        # aggregate intrusion detection.
                        if not intrusion:
                            # change the global matrix
                            self.glob[i][j][k] = 3
                            # change the local matrix
                            local_grids[li][lj][lk] = 3
                        else:
                            intru = True
                            break
                if intru:
                    break
            if intru:
                break

        # looping through all the elements in the local background grid of the current aggregate.
        # all elements in the blue dotted block.
        if not intru:
            for i in range(ib-1, Ib):
                for j in range(jb-1, Jb):
                    for k in range(kb-1, Kb):
                        # global element
                        global_element = i, j, k
                        # element central point (global background grid)
                        center_point = self._element_central_point(i, j, k, self.e)
                        # coordinates of the local background grid
                        li, lj, lk = (i-(ib-1), j-(jb-1), k-(kb-1))

                        if intru:
                            break

                        if local_grids[li][lj][lk] == 1:
                            #print("here")
                            if self._is_adjacent_to_aggregate(local_grids, li, lj, lk):
                                # These print statements are my little checkpoints.
                                #print(True)
                                # set the value of the local grid to 2 (representing ITZ)
                                local_grids[li][lj][lk] = 2
                                # set the corresponding value in the global grid to equal 2.
                                self.glob[i][j][k] = 2    

                    if intru:
                        break

                if intru:
                    break

        return local_grids, intru
    
        
# trying to generate and place 10 new aggregates.
n = Mesoscale(300, 300, 300, 4, 20, 20000)