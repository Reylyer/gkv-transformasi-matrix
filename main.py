from turtle import left
from manim import *
import numpy as np
from scipy.fftpack import shift

np.set_printoptions(threshold=np.inf)


class introNama(Scene):
    def construct(self):
        judul = Tex("Transfomasi")
        nama = [Tex("1. Faizal Husain Adiasha - 24060121140115"),
                Tex("2. Givandra Haikal Adjie - 24060121130063"),
                Tex("3. Zidan Rafindra Utomo - 24060121130051")
        ]

        nama[0].next_to(judul, DOWN).align_to(judul, OUT)
        nama[1].next_to(nama[0], DOWN).align_to(judul, OUT)
        nama[2].next_to(nama[1], DOWN).align_to(judul, OUT)

        tex_gr = VGroup(judul, *nama)
        tex_gr.move_to(ORIGIN)
        tex_gr.scale(1)
        judul.shift(UP)

        self.play(Write(judul))
        self.play(Write(nama[0]))
        self.play(Write(nama[1]))
        self.play(Write(nama[2]))
        self.wait()
        self.play(FadeOut(tex_gr))



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
        # self.do_scale(2)
        # self.wait()
        # self.do_scale(1/2)
        # self.wait()
        # self.do_rotation(90)
        # self.wait(2)
        # self.do_rotation(-90)
        self.wait(1)
        # self.do_abrrot([2, 1, 0], [3, 3, 1], 5)

        self.move_camera(phi=60*DEGREES, theta=-45*DEGREES, focal_distance=15, frame_center=self.main_obj.get_center(), run_time=2)



    def construct_axis(self):
        axis = ThreeDAxes()
        axis.set_z_index(3)
        labz = axis.get_z_axis_label(Tex("$y$"))
        laby = axis.get_y_axis_label(Tex("$z$"))
        labx = axis.get_x_axis_label(Tex("$x$"))

        self.play(FadeIn(axis), Write(labz), Write(laby), Write(labx))

    def setup_polyhedra(self):
        # som som major
        poly_points = [
            # x  z  y
            [ 3, 0, 0], # V0 kanan
            [ 0, 0, 2], # V1 atas
            [-3, 0, 0], # V2 kiri
            [ 0, 0,-2], # V3 bawah
            [ 0, 1, 0], # V4 keluar
            [ 0,-1, 0]  # V5 crot
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

class TranslasiScene(Scene):
    def construct(self):

        xyz_aksen = Matrix([
            [r"x'"],
            [r"y'"],
            [r"z'"],
            [1]
        ], element_alignment_corner=OUT)
        xyz = Matrix([
            [r"x"],
            [r"y"],
            [r"z"],
            [1]
        ], element_alignment_corner=OUT)
        step1_mat = Matrix([
            [r"x\cdot1 + y\cdot0 + z\cdot0 + T_{x}\cdot1"],
            [r"x\cdot0 + y\cdot1 + z\cdot0 + T_{y}\cdot1"],
            [r"x\cdot0 + y\cdot0 + z\cdot1 + T_{z}\cdot1"],
            [1]
        ], element_alignment_corner=OUT)
        step2_mat = Matrix([
            [r"x + T_{x}"],
            [r"y + T_{y}"],
            [r"z + T_{z}"],
            [1]
        ], element_alignment_corner=OUT)
        raw_transf_mat = Matrix([
            [ 1, 0, 0, r"T_{x}"],
            [ 0, 1, 0, r"T_{y}"],
            [ 0, 0, 1, r"T_{z}"],
            [ 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        tranf_mat = Matrix([
            [ 1, 0, 0, 2],
            [ 0, 1, 0, 3],
            [ 0, 0, 1, 2],
            [ 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        vertex_mat = Matrix([
            [ 3, 0, 0,-3, 0, 0],
            [ 0, 0, 2, 0, 0,-2],
            [ 0, 1, 0, 0,-1, 0],
            [ 1, 1, 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        vert_res_mat = Matrix([
            [ 5, 2, 2,-1, 2, 2],
            [ 3, 3, 5, 3, 3, 1],
            [ 2, 3, 2, 2, 1, 2],
            [ 1, 1, 1, 1, 1, 1],
        ], element_alignment_corner=OUT)

        up_eq_sign = Tex("=")
        down_eq_sign = Tex("=")
        judul = Tex("Translasi")
        label_tranf_mat = Tex("Matrix translasi").set_color(YELLOW)
        label_titik = Tex("Matrix titik").set_color(YELLOW)

        judul.to_corner(LEFT + UP)
        vertex_mat.next_to(tranf_mat, RIGHT)
        up_eq_sign.next_to(tranf_mat, LEFT)
        vert_res_mat.next_to(tranf_mat, DOWN)
        vert_res_mat.align_to(tranf_mat, LEFT)
        down_eq_sign.next_to(vert_res_mat, LEFT)
        xyz_aksen.next_to(up_eq_sign, LEFT)
        xyz.next_to(tranf_mat, RIGHT)
        step1_mat.next_to(tranf_mat, DOWN)
        step1_mat.align_to(tranf_mat, LEFT)
        step2_mat.next_to(tranf_mat, DOWN)
        step2_mat.align_to(tranf_mat, LEFT)

        label_tranf_mat.next_to(tranf_mat, UP)
        label_titik.next_to(vertex_mat, UP)


        all_group = VGroup(xyz_aksen, tranf_mat, vertex_mat, down_eq_sign, up_eq_sign, vert_res_mat, label_titik, label_tranf_mat, xyz, step1_mat, step2_mat)
        all_group.move_to(ORIGIN)
        all_group.scale(0.5)
        raw_transf_mat.scale(0.5)
        perkalian_group = VGroup(tranf_mat, vertex_mat)
        raw_transf_mat2 = raw_transf_mat.copy().move_to(tranf_mat)

        self.play(Write(judul))
        self.wait()
        self.play(Write(raw_transf_mat))
        self.wait()
        self.play(Transform(raw_transf_mat, tranf_mat))
        self.remove(tranf_mat)
        self.wait()
        self.play(Write(xyz_aksen), Write(up_eq_sign), Write(vertex_mat))
        self.wait()
        self.play(FadeIn(label_tranf_mat, shift=DOWN))
        self.wait
        self.play(FadeIn(label_titik, shift=DOWN))
        self.wait()
        self.play(FadeOut(label_tranf_mat, shift=UP), FadeOut(label_titik, shift=UP))
        self.wait()
        self.play(Write(down_eq_sign), TransformFromCopy(perkalian_group, vert_res_mat))
        self.wait()

        self.play(Transform(raw_transf_mat, raw_transf_mat2), vertex_mat.animate.become(xyz), FadeOut(vert_res_mat))
        self.wait()
        group_nyempil = VGroup(raw_transf_mat, vertex_mat)
        self.play(TransformFromCopy(group_nyempil, step1_mat))
        self.wait()
        self.play(Transform(step1_mat, step2_mat))
        self.wait()


class SkalaScene(Scene):
    def construct(self):

        xyz_aksen = Matrix([
            [r"x'"],
            [r"y'"],
            [r"z'"],
            [1]
        ], element_alignment_corner=OUT)
        xyz = Matrix([
            [r"x"],
            [r"y"],
            [r"z"],
            [1]
        ], element_alignment_corner=OUT)
        step1_mat = Matrix([
            [r"x \cdot S_{x} + y\cdot 0 + z\cdot 0 + 0 \cdot1"],
            [r"x \cdot 0 + y \cdot S_{y} + z\cdot 0 + 0 \cdot1"],
            [r"x \cdot 0 + y \cdot 0 + z \cdot S_{z} + 0 \cdot1"],
            [1]
        ], element_alignment_corner=OUT)
        step2_mat = Matrix([
            [r"x \cdot S_{x}"],
            [r"y \cdot S_{y}"],
            [r"z \cdot S_{z}"],
            [1]
        ], element_alignment_corner=OUT)
        raw_transf_mat = Matrix([
            [ r"S_{x}", 0, 0, 0],
            [ 0, r"S_{y}", 0, 0],
            [ 0, 0, r"S_{z}", 0],
            [ 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        tranf_mat = Matrix([
            [ 2, 0, 0, 0],
            [ 0, 2, 0, 0],
            [ 0, 0, 2, 0],
            [ 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        vertex_mat = Matrix([
            [ 3, 0, 0,-3, 0, 0],
            [ 0, 0, 2, 0, 0,-2],
            [ 0, 1, 0, 0,-1, 0],
            [ 1, 1, 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        vert_res_mat = Matrix([
            [ 6, 0, 0,-6, 0, 0],
            [ 0, 0, 4, 0, 0,-4],
            [ 0, 2, 0, 0,-2, 0],
            [ 1, 1, 1, 1, 1, 1],
        ], element_alignment_corner=OUT)

        up_eq_sign = Tex("=")
        down_eq_sign = Tex("=")
        judul = Tex("Skala")
        label_tranf_mat = Tex("Matrix skala").set_color(YELLOW)
        label_titik = Tex("Matrix titik").set_color(YELLOW)

        judul.to_corner(LEFT + UP)
        vertex_mat.next_to(tranf_mat, RIGHT)
        up_eq_sign.next_to(tranf_mat, LEFT)
        vert_res_mat.next_to(tranf_mat, DOWN)
        vert_res_mat.align_to(tranf_mat, LEFT)
        down_eq_sign.next_to(vert_res_mat, LEFT)
        xyz_aksen.next_to(up_eq_sign, LEFT)
        xyz.next_to(tranf_mat, RIGHT)
        step1_mat.next_to(tranf_mat, DOWN)
        step1_mat.align_to(tranf_mat, LEFT)
        step2_mat.next_to(tranf_mat, DOWN)
        step2_mat.align_to(tranf_mat, LEFT)

        label_tranf_mat.next_to(tranf_mat, UP)
        label_titik.next_to(vertex_mat, UP)


        all_group = VGroup(xyz_aksen, tranf_mat, vertex_mat, down_eq_sign, up_eq_sign, vert_res_mat, label_titik, label_tranf_mat, xyz, step1_mat, step2_mat)
        all_group.move_to(ORIGIN)
        all_group.scale(0.5)
        raw_transf_mat.scale(0.5)
        perkalian_group = VGroup(tranf_mat, vertex_mat)
        raw_transf_mat2 = raw_transf_mat.copy().move_to(tranf_mat)

        self.play(Write(judul))
        self.wait()
        self.play(Write(raw_transf_mat))
        self.wait()
        self.play(Transform(raw_transf_mat, tranf_mat))
        self.remove(tranf_mat)
        self.wait()
        self.play(Write(xyz_aksen), Write(up_eq_sign), Write(vertex_mat))
        self.wait()
        self.play(FadeIn(label_tranf_mat, shift=DOWN))
        self.wait
        self.play(FadeIn(label_titik, shift=DOWN))
        self.wait()
        self.play(FadeOut(label_tranf_mat, shift=UP), FadeOut(label_titik, shift=UP))
        self.wait()
        self.play(Write(down_eq_sign), TransformFromCopy(perkalian_group, vert_res_mat))
        self.wait()

        self.play(Transform(raw_transf_mat, raw_transf_mat2), vertex_mat.animate.become(xyz), FadeOut(vert_res_mat))
        self.wait()
        group_nyempil = VGroup(raw_transf_mat, vertex_mat)
        self.play(TransformFromCopy(group_nyempil, step1_mat))
        self.wait()
        self.play(Transform(step1_mat, step2_mat))
        self.wait()


class RotasiScene(Scene):
    def construct(self):
        xyz_aksen = Matrix([
            [r"x'"],
            [r"y'"],
            [r"z'"],
            [1]
        ], element_alignment_corner=OUT)
        xyz = Matrix([
            [r"x"],
            [r"y"],
            [r"z"],
            [1]
        ], element_alignment_corner=OUT)
        step1_mat = Matrix([
            [r"x\cdot1 + y\cdot0 + z\cdot0 + T_{x}\cdot1"],
            [r"x\cdot0 + y\cdot1 + z\cdot0 + T_{y}\cdot1"],
            [r"x\cdot0 + y\cdot0 + z\cdot1 + T_{z}\cdot1"],
            [1]
        ], element_alignment_corner=OUT)
        step2_mat = Matrix([
            [r"x + T_{x}"],
            [r"y + T_{y}"],
            [r"z + T_{z}"],
            [1]
        ], element_alignment_corner=OUT)
        raw_transf_mat = Matrix([
            [ 1, 0, 0, r"T_{x}"],
            [ 0, 1, 0, r"T_{y}"],
            [ 0, 0, 1, r"T_{z}"],
            [ 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        tranf_mat = Matrix([
            [ 1, 0, 0, 2],
            [ 0, 1, 0, 3],
            [ 0, 0, 1, 3],
            [ 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        vertex_mat = Matrix([
            [ 3, 0, 0,-3, 0, 0],
            [ 0, 0, 2, 0, 0,-2],
            [ 0, 1, 0, 0,-1, 0],
            [ 1, 1, 1, 1, 1, 1],
        ], element_alignment_corner=OUT)
        vert_res_mat = Matrix([
            [ 5, 2, 2,-1, 2, 2],
            [ 3, 3, 5, 3, 3, 1],
            [ 2, 3, 2, 2, 1, 2],
            [ 1, 1, 1, 1, 1, 1],
        ], element_alignment_corner=OUT)

        up_eq_sign = Tex("=")
        down_eq_sign = Tex("=")
        judul = Tex("Translasi")
        label_tranf_mat = Tex("Matrix transformasi").set_color(YELLOW)
        label_titik = Tex("Matrix titik").set_color(YELLOW)

        judul.to_corner(LEFT + UP)
        vertex_mat.next_to(tranf_mat, RIGHT)
        up_eq_sign.next_to(tranf_mat, LEFT)
        vert_res_mat.next_to(tranf_mat, DOWN)
        vert_res_mat.align_to(tranf_mat, LEFT)
        down_eq_sign.next_to(vert_res_mat, LEFT)
        xyz_aksen.next_to(up_eq_sign, LEFT)
        xyz.next_to(tranf_mat, RIGHT)
        step1_mat.next_to(tranf_mat, DOWN)
        step1_mat.align_to(tranf_mat, LEFT)
        step2_mat.next_to(tranf_mat, DOWN)
        step2_mat.align_to(tranf_mat, LEFT)

        label_tranf_mat.next_to(tranf_mat, UP)
        label_titik.next_to(vertex_mat, UP)


        all_group = VGroup(xyz_aksen, tranf_mat, vertex_mat, down_eq_sign, up_eq_sign, vert_res_mat, label_titik, label_tranf_mat, xyz, step1_mat, step2_mat)
        all_group.move_to(ORIGIN)
        all_group.scale(0.5)
        raw_transf_mat.scale(0.5)
        perkalian_group = VGroup(tranf_mat, vertex_mat)
        raw_transf_mat2 = raw_transf_mat.copy().move_to(tranf_mat)

        self.play(Write(judul))
        self.wait()
        self.play(Write(raw_transf_mat))
        self.wait()
        self.play(Transform(raw_transf_mat, tranf_mat))
        self.remove(tranf_mat)
        self.wait()
        self.play(Write(xyz_aksen), Write(up_eq_sign), Write(vertex_mat))
        self.wait()
        self.play(FadeIn(label_tranf_mat, shift=DOWN))
        self.wait
        self.play(FadeIn(label_titik, shift=DOWN))
        self.wait()
        self.play(FadeOut(label_tranf_mat, shift=UP), FadeOut(label_titik, shift=UP))
        self.wait()
        self.play(Write(down_eq_sign), TransformFromCopy(perkalian_group, vert_res_mat))
        self.wait()

        self.play(Transform(raw_transf_mat, raw_transf_mat2), vertex_mat.animate.become(xyz), FadeOut(vert_res_mat))
        self.wait()
        group_nyempil = VGroup(raw_transf_mat, vertex_mat)
        self.play(TransformFromCopy(group_nyempil, step1_mat))
        self.wait()
        self.play(Transform(step1_mat, step2_mat))
        self.wait()


class TestScene(Scene):
    def construct(self):
        r = Rectangle()
        r2 = Rectangle(height=4, width=1)
        r2.shift(3*LEFT)
        r.shift(UP)
        self.play(r.animate.shift(RIGHT*2), rate_func=linear)
        self.play(r.animate.shift(DOWN*2), rate_func=smooth)
        self.play(Transform(r, r2))
        self.wait()
