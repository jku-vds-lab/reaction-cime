import math
import time

import numba
import numpy as np
from sklearn.neighbors import KDTree


@numba.jit(nopython=True, parallel=False)
def _density_calculation(count_map, mask, mask_size, x_min, x_max, y_min, y_max, nr_columns, nr_rows):
    dummy_points_candidates = []

    for row in range(nr_rows):
        y_ = row * (y_max - y_min) / (nr_rows - 1) + y_min

        for column in range(nr_columns):
            if count_map[row][column] == 0:
                x_ = column * (x_max - x_min) / (nr_columns - 1) + x_min
                density = 0

                for i in range(mask_size):
                    for j in range(mask_size):
                        r = row - (int(mask_size / 2)) + j
                        c = column - (int(mask_size / 2)) + i

                        if (0 <= r < nr_rows) and (0 <= c < nr_columns):
                            density += mask[i][j] * count_map[r][c]

                dummy_points_candidates.append([x_, y_, density, -1])

    return dummy_points_candidates


class DGrid:
    def __init__(self, icon_width=1, icon_height=1, delta=None, callbacks=[]):
        self.icon_width_ = icon_width
        self.icon_height_ = icon_height
        self.delta_ = delta

        if self.delta_ is None:
            self.delta_ = 1

        self.grid_ = []
        self.callbacks = callbacks

    def _trigger_callbacks(self, msg):
        if self.callbacks is not None:
            for callback in self.callbacks:
                callback(msg)

    def _fit(self, y):
        self._trigger_callbacks("start fit")
        # calculating the bounding box
        max_coordinates = np.amax(y, axis=0)
        min_coordinates = np.amin(y, axis=0)
        bounding_box_width = max_coordinates[0] - min_coordinates[0]
        bounding_box_height = max_coordinates[1] - min_coordinates[1]

        # defining the number of rows and columns
        nr_columns = math.ceil((self.delta_ * bounding_box_width) / self.icon_width_)
        nr_rows = math.ceil((self.delta_ * bounding_box_height) / self.icon_height_)

        # if the number of rows and columns are not enough to fit all data instances, increase delta
        if nr_rows * nr_columns < len(y):
            nr_columns = bounding_box_width / self.icon_width_
            nr_rows = bounding_box_height / self.icon_height_
            self.delta_ = math.sqrt(len(y) / (nr_rows * nr_columns))
            nr_columns = math.ceil(self.delta_ * nr_columns)
            nr_rows = math.ceil(self.delta_ * nr_rows)

            print(
                "There is not enough space to remove overlaps! Setting delta to {0}, the smallest possible number "
                "to fully remove overlaps. Increase it if more empty space is required.".format(self.delta_)
            )

        # add the original points
        def to_grid_cell(id_, x_, y_):
            return {"id": id_, "x": x_, "y": y_, "i": 0, "j": 0, "dummy": False}

        self._trigger_callbacks("refactor data")
        for i in range(len(y)):
            self.grid_.append(to_grid_cell(i, y[i][0], y[i][1]))

        # add the dummy points
        start_time = time.time()
        self._add_dummy_points(min_coordinates[0], max_coordinates[0], min_coordinates[1], max_coordinates[1], nr_columns, nr_rows)
        print("--- Add dummy points executed in %s seconds ---" % (time.time() - start_time))

        # execute
        self._trigger_callbacks("grid assignment")
        start_time = time.time()
        self.grid_ = DGrid._grid_rec(self.grid_, nr_rows, nr_columns, 0, 0)
        self.grid_.sort(key=lambda v: v.get("id"))
        print("--- Grid assignment executed in %s seconds ---" % (time.time() - start_time))

        self._trigger_callbacks("transform")
        # returning the overlap free scatterplot
        transformed = []
        for i in range(len(self.grid_)):
            if self.grid_[i]["dummy"] is False:
                transformed.append(
                    np.array(
                        [
                            min_coordinates[0] + self.grid_[i]["j"] * self.icon_width_,
                            min_coordinates[1] + self.grid_[i]["i"] * self.icon_height_,
                        ]
                    )
                )

        self._trigger_callbacks("finish")
        return np.array(transformed)

    def fit_transform(self, y):
        return self._fit(y)

    def fit(self, y):
        return self._fit(y)

    @staticmethod
    def _split_grid(grid, cut_point, direction):
        if direction == "x":
            grid.sort(key=lambda cel: (cel["x"], cel["y"]))
        else:
            grid.sort(key=lambda cel: (cel["y"], cel["x"]))

        grid0 = grid[:cut_point]
        grid1 = []
        if cut_point < len(grid):
            grid1 = grid[-(len(grid) - cut_point) :]

        return grid0, grid1

    @staticmethod
    def _grid_rec(grid, r, s, i, j):
        size = len(grid)

        if size > 0:
            if size == 1:
                grid[0]["i"] = i
                grid[0]["j"] = j
            else:
                if r > s:
                    half_rows = int(math.ceil(r / 2.0))
                    grid0, grid1 = DGrid._split_grid(grid, min(size, half_rows * s), "y")
                    DGrid._grid_rec(grid0, half_rows, s, i, j)
                    DGrid._grid_rec(grid1, (r - half_rows), s, (i + half_rows), j)
                else:
                    half_columns = int(math.ceil(s / 2.0))
                    grid0, grid1 = DGrid._split_grid(grid, min(size, half_columns * r), "x")
                    DGrid._grid_rec(grid0, r, half_columns, i, j)
                    DGrid._grid_rec(grid1, r, (s - half_columns), i, (j + half_columns))

        return grid

    def _add_dummy_points(self, x_min, x_max, y_min, y_max, nr_columns, nr_rows):
        self._trigger_callbacks("add dummy points")
        size = len(self.grid_)

        # counting the number of points per grid cell
        count_map = np.zeros((nr_rows, nr_columns), dtype=np.uint32)

        for i in range(size):
            col = math.ceil(((self.grid_[i]["x"] - x_min) / (x_max - x_min)) * (nr_columns - 1))
            row = math.ceil(((self.grid_[i]["y"] - y_min) / (y_max - y_min)) * (nr_rows - 1))
            count_map[row][col] = count_map[row][col] + 1

        # calculating the gaussian mask
        mask_size = int(max(3, ((x_max - x_min) * (y_max - y_min)) / (size * self.icon_width_ * self.icon_height_)))
        mask_size = mask_size + 1 if mask_size % 2 == 0 else mask_size
        mask = DGrid._gaussian_mask(mask_size, (mask_size - 1) / 6.0)

        self._trigger_callbacks("density calculation")
        # creating all dummy candidates
        dummy_points_candidates = _density_calculation(count_map, mask, mask_size, x_min, x_max, y_min, y_max, nr_columns, nr_rows)

        self._trigger_callbacks("post process dummy points")
        # sorting candidates using density
        dummy_points_candidates.sort(key=lambda x: x[2])

        # defining the number of required dummy points
        nr_dummy_points = min((nr_rows * nr_columns) - size, len(dummy_points_candidates))

        # checking if density is not enough to decide the correct dummy points
        if (
            len(dummy_points_candidates) > nr_dummy_points
            and math.fabs(dummy_points_candidates[nr_dummy_points - 1][2] - dummy_points_candidates[nr_dummy_points][2]) < 0.0001
        ):

            # if not, create a kd-tree to find the nearest point in the original layout
            original_points = []
            for i in range(size):
                # adding the original points
                x_ = self.grid_[i]["x"]
                y_ = self.grid_[i]["y"]
                original_points.append([x_, y_])

            tree = KDTree(original_points, leaf_size=2)

            # add the distance information for the "undecided" dummy points
            for i in range(len(dummy_points_candidates)):
                if math.fabs(dummy_points_candidates[nr_dummy_points - 1][2] - dummy_points_candidates[i][2]) < 0.0001:
                    dummy_points_candidates[i][3] = float(
                        tree.query([[dummy_points_candidates[i][0], dummy_points_candidates[i][1]]], 1)[0]
                    )

            # sort the candidates again using density and distance
            dummy_points_candidates.sort(key=lambda x: (x[2], x[3]))

        for i in range(nr_dummy_points):
            self.grid_.append(
                {"id": size + i, "x": dummy_points_candidates[i][0], "y": dummy_points_candidates[i][1], "i": 0, "j": 0, "dummy": True}
            )

        return

    @staticmethod
    def _gaussian_mask(size, sigma):
        mask = np.zeros((size, size), dtype=np.float32)

        for i in range(size):
            y = int(i - int(size / 2))

            for j in range(size):
                x = int(j - int(size / 2))
                mask[i][j] = 1.0 / (2 * math.pi * sigma * sigma) * math.exp(-(x * x + y * y) / (2 * sigma * sigma))

        return mask
