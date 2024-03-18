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