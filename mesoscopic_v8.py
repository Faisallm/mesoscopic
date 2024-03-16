import numpy as np
import random
import math

class Mesoscale:
    
    
    def __init__(self, L, W, H, dmin, dmax, number):
        
        
        self.L = L
        self.W = W
        self.H = H
        self.dmin = dmin
        self.dmax = dmax
        self.number = number
        
        # mesh size
        self.e = (1/5)*self.dmin
        
        # number of elements along x,  y, z
        self.l = int(self.L//self.e)
        self.m = int(self.W//self.e)
        self.n = int(self.H//self.e)
        self.total_elements = self.l * self.m * self.n
        
        # test to ensure dimensions and mesh size are consistent
        assert((self.L*self.W*self.H) == (self.l*self.e * self.m * self.e * self.n * self.e))
        
        # material attribute for global matrix
        self.glob = self._global_matrix(self.l, self.m, self.n)
        
        # aggregate size ranges
        self.aggregate_size_ranges = [[5, 10], [10, 15], [15, 20]]
#         self.aggregate_size_ranges = [[0.15, 0.300], [0.300, 0.6], [0.6, 1.18],
#             [1.18, 2.36], [2.36, 4.75], [4.75, 9.5]]
#         self.aggregate_size_ranges = [[4.75, 9.5], [9.5, 12.5], [12.5, 20.0]]
    
        self.volume_fractions = [self._aggregate_volume_fraction(i) for i in self.aggregate_size_ranges]
#         print(self.volume_fractions)
        
        self.agg_fraction = [int(i * self.number) for i in self.volume_fractions]
#         print(self.agg_fraction)
        
        self.aggs = []
        
        for i, size in enumerate(self.aggregate_size_ranges):
            # we are generating a certain num of aggregates for each size_range
            aggu = [self._generate_polyhedron(size[0], size[1]) for _ in range(self.agg_fraction[i])]
            
            for agger in aggu:
                self.aggs.append(agger)
                
            
#         print(f"Number of generated aggregates: {len(self.aggs)}")
#         print(self.aggs)
        
#         self.probabilities = [self._fuller_curve(i) for i in self.aggregate_size_ranges]
        
#         self.aggs = [self._generate_polyhedron(dmin, dmax) \
#                      for dmin, dmax in random.choices([i for i in self.aggregate_size_ranges], weights=self.probabilities) \
#                      for j in range(self.number)]
        
        # aggregates arranged in order of decreasing size.
        self.arranged_aggs = sorted(self.aggs, key=lambda x: self._polyhedral_area(x), reverse=True)
        
#         agg1 = self.arranged_aggs[0]
        
#         self.translated_agg = self._translate(agg1,
#                                           self._random_translation_point(self.l, self.m, self.n, self.e))
        
        # testing
        self.agg_counts = 0
        for i, a in enumerate(self.arranged_aggs):
            self._condition = True
            self._count = 0
            while self._condition and self._count < 20:
                self.agg = self._translate(a,
                                          self._random_translation_point(self.l, self.m, self.n, self.e))
                # local grid
                # If their is an intrusion, I want the aggregate to be randomly translated to a new
                # location for placing for atleast 500 times.
                self.loc, self._condition = self._identify_aggregate_and_itz(self.agg)
                print(f"Local: {self.loc}")
                # checking if the aggregate was placed
                if not self._condition:
                    # break out of loop
                    self._condition = False
                    self._count = 0
                    self.agg_counts += 1
                    print("no intrusion")
                else:
                    self._count += 1;
                    print("intrusion detected!")
            # to carry out the same check for other aggregates
            
            if not self._condition:
                print(f"Placed #{self.agg_counts} aggregate.")
                
            self._condition = True
            
        
        
        #print(f"Local: {self.loc}")
        print(f"Global: {self.glob}")
        
    
    def _aggregate_volume_fraction(self, aggregate_size_range):
        """This function calculates the volume fraction of a 
        particular aggregate size range (di+1, di).
        
        The aggregate content of concrete is closely related to the size
        distribution of aggregates, which affect the macroscopic mechanical
        properties of concrete."""
        
        size1, size2 = aggregate_size_range[0], aggregate_size_range[1]
        
        pdii, pdi = self._fuller(size1), self._fuller(size2)
        
        pdmax, pdmin = self._fuller(self.dmax), self._fuller(self.dmin)
#         print(((pdii - pdi) / (pdmax - pdmin))*-1)
        
        return ((pdii - pdi) / (pdmax - pdmin))*-1
    
    def _fuller(self, d):
        # these values do not change.
        n = 0.5
        # upper limit of aggregate size
        D = 20
        return ((d/D)**n)
        
        
    def _minimum_and_maximum(self, agg):
        """According to the new coordinates of the vertices of the...
        polyhedral aggregate. The minimum values (Xmin, Ymin, Zmin) and
        the maximum values (Xmax, Ymax, Zmax) of all the
        vertex coordinates in the three-dimensional direction
        are obtained in this function."""

        # axis 0 refers to the column.
        # minimum and maxmimum values of the vertices of an aggregate.
        x_min, y_min, z_min = np.amin(agg, axis=0)
        x_max, y_max, z_max = np.amax(agg, axis=0)

        return (int(x_min), int(y_min), int(z_min)), (int(x_max), int(y_max), int(z_max))
    
    
    def _element_central_point(self, i: int, j: int, k: int, e):
        """The coordinates of the central point of any element
        (i, j, k) in the global background grid.
        The type of array returned is a numpy array to optimize
           numerical calculations."""
        
        return np.array(((i*e) + e/2, (j*e)+e/2, (k*e)+e/2))
    
    
    def _section_global_background_grid(self, aggregate):
        """The bounding box ((ir, jr, kr), (Ir, Jr, Kr)) of the newly placed aggregate can be
        determined using(red dashes) without considering the ITZ.

        While the coordinates ((ib, jb, kb), (Ib, Jb, Kb)) represents the 
        bounding box without considering the ITZ."""


        ((x_min, y_min, z_min), (x_max, y_max, z_max)) = self._minimum_and_maximum(aggregate)

        # bounding box of the newly placed aggregate (consider aggregate elemennts).
        (ir, jr, kr) = np.max([0, x_min//self.e]), np.max([0, y_min//self.e]), np.max([0, z_min//self.e])
        (Ir, Jr, Kr) = np.min([self.l, x_max//self.e]), np.min([self.m, y_max//self.e]), np.min([self.n, z_max//self.e])

        # bounding box of the newly placed aggregate considering ITZ elements
        (ib, jb, kb) = np.max([0, (x_min//self.e)-1]), \
            np.max([0, (y_min//self.e)-1]), np.max([0, (z_min//self.e)-1])
        
        (Ib, Jb, Kb) = np.min([self.l, (x_max//self.e)+1]),\
            np.min([self.m, (y_max//self.e)+1]), np.min([self.n, (z_max//self.e)+1])

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
    
    
    # I now agree that this function is correct.
    # it makes sense now.
    def _local_grid(self, lower_limit, upper_limit):
        """
        We have to generate the local background grid from the global
        background grid.

        Always remember that range excludes the outer limit, so don't
        forget to add + 1.
        All elements in the local background grid are initialized
        as mortar elements.
        """
        # they are being mapped to int.
        (ib, jb, kb), (Ib, Jb, Kb) = map(int, lower_limit), map(int, upper_limit)

        ooo_array = []
        for i in range(ib, Ib+1):
            oo_array = []
            for j in range(jb, Jb+1):
                o_array = []
                for k in range(kb, Kb+1):
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
    
    # another area of concern!!!
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
        ir, jr, kr = int(ir), int(jr), int(kr)
        Ir, Jr, Kr = int(Ir), int(Jr), int(Kr)
        
        # this function uses alot of compute unneccessarily (hence the introduction of breaks).
        for i in range(ir-1, Ir):
            for j in range(jr-1, Jr):
                for k in range(kr-1, Kr):
                    if self.glob[i][j][k] == 3 or self.glob[i][j][k] == 2: 
                        contains = True
                        break
                if contains:
                    break
            if contains:
                break
                        
        return contains
    
    
    def _identify_aggregate_and_itz(self, agg):
        """
        1, loop through all the elements in the local background grid of current aggregate.
        i.e all the elements in the red dotted block.

        2, According to the element number of the local background grid,
        the corresponding element of the global background grid is obtained.
        and then the coordinates of the center points are obtained (eq 4).
        """

        # a section of the global matrix
        ((ir, jr, kr), (Ir, Jr, Kr)), ((ib, jb, kb), (Ib, Jb, Kb)) = self._section_global_background_grid(agg)
        local_grids = self._local_grid((ib, jb, kb), (Ib, Jb, Kb))
#         print(local_grids.shape)
        intrusion = self._aggregate_intrusion(ib, Ib, jb, Jb, kb, Kb)
        intru = False

#         ir, jr, kr = int(ir), int(jr), int(kr)
#         Ir, Jr, Kr = int(Ir), int(Jr), int(Kr)
        ib, jb, kb = int(ib), int(jb), int(kb)
        Ib, Jb, Kb = int(Ib), int(Jb), int(Kb)
        for i in range(ib, Ib):
            for j in range(jb, Jb):
                for k in range(kb, Kb):
                    global_element = i, j, k
                    center_point = self._element_central_point(i, j, k, self.e)
                    li, lj, lk = (i-ib, j-jb, k-kb)

                    # if the element is aggregate
                    if self._point_in_polyhedron(center_point, agg):
                        
                        # aggregate intrusion detection.
                        if not intrusion:
                            # change the global matrix
                            self.glob[i][j][k] = 3
                            # change the local matrix
                            # how do I know that its changing the correct position in the global matrix?
                            local_grids[li][lj][lk] = 3 
                        else:
                            intru = True
                            break
                            
#                     elif local_grids[li][lj][lk] == 1 and self._is_adjacent_to_aggregate(local_grids, li, lj, lk):
#                         # change element in local to ITZ.
#                         local_grids[li][lj][lk] = 2
#                         # change corresponding element in global to ITZ.
#                         self.glob[i][j][k] = 2
#                     else:
#                         # change element in local to ITZ.
#                         local_grids[li][lj][lk] = 1
#                         # change corresponding element in global to ITZ.
#                         self.glob[i][j][k] = 1
                        
                if intru:
                    break
            if intru:
                break

#         ib, jb, kb = int(ib), int(jb), int(kb)
#         Ib, Jb, Kb = int(Ib), int(Jb), int(Kb)
        if not intru:
            for i in range(ib, Ib):
                for j in range(jb, Jb):
                    for k in range(kb, Kb):
                        # coordinates of the global grid
                        global_element = i, j, k
                        # coordinates of the local background grid
                        li, lj, lk = (i-ib, j-jb, k-kb)
#                         print("itz")
#                         print(f"global: i: {i}, j:{j}, k: {k}")
#                         print(f"local: li: {li}, lj:{lj}, lk: {lk}")

                        if intru:
                            break

                        if local_grids[li][lj][lk] == 1:
                            if self._is_adjacent_to_aggregate(local_grids, li, lj, lk):
                                # change element in local to ITZ.
                                local_grids[li][lj][lk] = 2
                                # change corresponding element in global to ITZ.
                                self.glob[i][j][k] = 2   
                    if intru:
                        break

                if intru:
                    break

        return local_grids, intru
    
    
    def _locate_element(self, x: int, y: int, z: int, e):
        """If the coordinates of any point in the background
        mesh block are (x, y, z), the number (i, j, k) of the element
        where the point can be located is defined in this function.
        The type of array returned is a numpy array to optimize
           numerical calculations."""
        
        return np.array((x // e, y // e, z // e))
    
        
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
        
        n7 = random.random()
        n8 = random.random()
        n9 = random.random()

        xp, yp, zp = (n7 * l * e, n8 * m * e, n9 * n * e)

        return np.array((xp, yp, zp))
    
    
    def _translate(self, aggregate, translation_point):
        """The newly generated random polyhedral aggregate is directly put into
        the background grid space by RANDOM TRANSLATION.
        Assuming the coordinates of the random translation point are
        (xp, yp, zp)=(n7le, n8me, n9ne).
        Therefore, the coordinates of the vertices after translation are:
        Xi = xi + xp, Yi = yi + yp, Zi = zi + zp."""
        
        condition = True
        
        while condition:

            # generate a random translation point.
            xp, yp, zp = self._random_translation_point(self.l, self.m, self.n, self.e)
            # translate the aggregate to the new point.
            translated_aggregate = [(i[0] + xp, i[1] + yp, i[2] + zp) for i in aggregate]
            
            for tran_agg in translated_aggregate:
                
                element = self._locate_element(tran_agg[0], tran_agg[1], tran_agg[2], self.e)
                i, j, k = list(element)
                
                if not ((0 <= tran_agg[0] <= self.l * self.e) and 
                        (0 <= tran_agg[1] <= self.n * self.e) and 
                        (0 <= tran_agg[2] <= self.m * self.e)):
                    # if one of the coordinate is not located within the concrete specimen.
                    # break out of the loop, no need to continue
                    condition = True
                    break
                    
                elif self.glob[int(i)][int(j)][int(k)] != 1:
                    # even if the random translation point for the...
                    # aggregate is within the concrete specimen
                    # but located in an ITZ or aggregate element.
                    # generate a new translation point.
                    # break out of loop no need to continue
                    condition = True
                    break
                    
                else:
                    # here we will not break out
                    # condition has to be false for all the coordinates...
                    # of the translated aggregate for the translated...
                    # aggregate to be selected.
                    condition = False
                    
                    
            # get the element where this point is located
            #element = self._locate_element(xp, yp, zp, self.e)
            
            #i, j, k = list(element)
            
            # ensure random point P is mortar with (1)
            #if not self.glob[int(i)][int(j)][int(k)] == 1:
            #    condition = True

            
        return np.array(translated_aggregate)
    
        
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
    
        
    def _fuller_curve(self, d):
        # these values do not change.
        d=d[1]
        n = 0.5
        # upper limit of aggregate size
        D = 20
        return ((d/D)**n)
        
        
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

        The meso-scale finite element of concrete will be generated using this matrix.
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
        
        