from manim import *
import numpy as np

np.set_printoptions(threshold=np.inf)

class Transformasiabrrot(ThreeDScene, MovingCamera):
    def __init__(self):
        super().__init__()

    def construct(self):
        self.set_camera_orientation(phi=60*DEGREES, theta=-45*DEGREES, focal_distance=15)

        self.construct_axis()
        self.setup_polyhedra()
        self.write_mat()
        self.do_translate(2, 2, 3)
        self.wait()
        self.do_translate(-2, -2, -3)
        self.do_scale(2)
        self.wait()
        self.do_scale(1/2)
        self.wait()
        self.do_rotation(90)
        self.wait(2)
        self.do_rotation(-90)
        self.wait(1)
        # self.do_abrrot([2, 1, 0], [3, 3, 1], 5)

        self.move_camera(phi=60*DEGREES, theta=-45*DEGREES, focal_distance=15, frame_center=self.main_obj.get_center(), run_time=2)



    def construct_axis(self):
        axis = ThreeDAxes()
        axis.set_z_index(3)
        labz = axis.get_z_axis_label(Tex("$z$"))
        laby = axis.get_y_axis_label(Tex("$y$"))
        labx = axis.get_x_axis_label(Tex("$x$"))

        self.play(FadeIn(axis), Write(labz), Write(laby), Write(labx))

    def setup_polyhedra(self):
        # som som major
        poly_points = [
            [ 3, 0, 0], # V0 kanan
            [ 0, 2, 0], # V1 atas
            [-3, 0, 0], # V2 kiri
            [ 0,-2, 0], # V3 bawah
            [ 0, 0, 1], # V4 keluar
            [ 0, 0,-1]  # V5 crot
        ]
        faces_list = [
            [0, 1, 4],
            [0, 1, 5],
            [1, 2, 4],
            [1, 2, 5],
            [2, 3, 4],
            [2, 3, 5],
            [3, 0, 4],
            [3, 0, 5],
        ]
        self.main_obj = main_obj = Polyhedron(vertex_coords=poly_points, faces_list=faces_list)
        self.play(DrawBorderThenFill(main_obj), run_time=2)
        self.move_camera(phi=60*DEGREES, theta=-45*DEGREES, focal_distance=15, frame_center=main_obj.get_center(), run_time=2)
        self.wait(0.5)

    def write_mat(self):

        def get_poly_coords():
            coords_faces = self.main_obj.extract_face_coords()
            # [print(n, n.shape) for n in coords_faces[0] + coords_faces[5]]
            coords = np.stack(coords_faces[0] + coords_faces[5])
            coords = coords.round(2)
            # print(coords)
            coords = coords.transpose()
            coords = np.append(coords, [[1, 1, 1, 1, 1, 1]], axis=0)
            # print(coords)
            return coords

        def matrix_updater(mob:Matrix):
            # https://www.reddit.com/r/manim/comments/oid6hv/comment/h4vxk5y/?utm_source=share&utm_medium=web2x&context=3
            newMat = Matrix(get_poly_coords(), h_buff=2)
            newMat.scale(0.5)
            newMat.to_corner(UP + LEFT)
            newMat.shift(OUT*20)
            mob.become(newMat, copy_submobjects=False)

        pmat = Matrix(get_poly_coords(), h_buff=2)
        self.add_fixed_in_frame_mobjects(pmat)
        pmat.scale(0.5)
        pmat.to_corner(UP + LEFT)

        self.play(Write(pmat))
        # self.add(pmat)
        self.wait()
        pmat.add_updater(matrix_updater)
    
    def do_abrrot(self, startp, endp, dur):

        def main_obj_rot_updater(main_obj_r, dt):
            main_obj_r.rotate(angle=dt, axis=AB, about_point=point_rot)

        # https://gamedev.stackexchange.com/questions/72528/how-can-i-project-a-3d-point-onto-a-3d-line
        self.rot_axis = rot_axis = Line3D(startp, endp)
        AB = rot_axis.end - rot_axis.start
        AP = self.main_obj.get_center() - rot_axis.start
        point_rot = rot_axis.start + np.dot(AP, AB) / np.dot(AB, AB) * AB

        self.play(Write(rot_axis))
        # self.add(rot_axis)

        self.move_camera(phi=60*DEGREES, theta=-45*DEGREES, focal_distance=20, frame_center=rot_axis.get_center())
        self.begin_ambient_camera_rotation(45*DEGREES/3, about='theta')
        self.main_obj.add_updater(main_obj_rot_updater)
        self.wait(dur)
        self.main_obj.remove_updater(main_obj_rot_updater)
        self.stop_ambient_camera_rotation(about="theta")
    
    def do_translate(self, x, y, z):
        self.play(self.main_obj.animate.shift(x*RIGHT + y*UP + z*OUT), run_time=2)
        self.move_camera(phi=60*DEGREES, theta=-45*DEGREES, focal_distance=15, frame_center=self.main_obj.get_center(), run_time=2)

    def do_scale(self, scale_fact):
        phantom_axes = ThreeDAxes()
        phantom_axes.set_opacity(0.0)
        scl_group = VGroup(phantom_axes, self.main_obj)
        self.play(scl_group.animate.scale(scale_fact), run_time=2)
        self.move_camera(phi=60*DEGREES, theta=-45*DEGREES, focal_distance=15, frame_center=self.main_obj.get_center(), run_time=2)

    def do_rotation(self, deg):
        self.play(self.main_obj.animate.rotate(angle=deg*DEGREES, axis=UP, about_point=ORIGIN), run_time=2)


class introNama(Scene):
    def construct(self):
        judul = Tex("Transfomasi")
        nama = [Tex("Faizal Husain Adiasha"),
                Tex("Givandra Haikal Adjie"),
                Tex("Zidan Rafindra Utomo")
        ]

        judul.to_edge(UP)
        judul.shift(LEFT)
        nama[0].next_to(judul, DOWN).align_to(judul, OUT)
        nama[1].next_to(nama[0], DOWN).align_to(judul, OUT)
        nama[2].next_to(nama[1], DOWN).align_to(judul, OUT)

        tex_gr = VGroup(judul, *nama)
        tex_gr.move_to(ORIGIN)
        tex_gr.scale(1.5)

        self.play(Write(judul))
        self.play(Write(nama[0]))
        self.play(Write(nama[1]))
        self.play(Write(nama[2]))
        self.wait()
        self.play(FadeOut(tex_gr))
