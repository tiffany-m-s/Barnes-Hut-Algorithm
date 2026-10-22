import os, copy
from dataclasses import dataclass
from datatypes import OrderedPair, Star, Universe, Quadrant, Node, QuadTree, center_of_gravity
from engine import generate_quadtree

@dataclass
class Center_Of_Gravity_Test:
    stars: list[Star]
    result: OrderedPair

@dataclass
class Create_Children_Test:
    node: Node
    result: list[Node]

@dataclass
class Find_Child_Test:
    star: Star
    node: Node
    result: Node

@dataclass
class Node_Insert_Test:
    root: Node
    star: Star
    result: Node

@dataclass
class Generate_QuadTree_Test:
    universe: Universe
    result: Node


def read_center_of_gravity_sample(file):
    stars = []
    with open(file, "r") as f:
        for line in f:
            # The each line is a star's attributes
            star = Star()
            star_attributes = line.rstrip("\n").split(" ")
            star.position = OrderedPair(float(star_attributes[0]), float(star_attributes[1]))
            star.mass = float(star_attributes[2])
            stars.append(star)
    return stars

def read_ordered_pair_from_file(file):
    with open(file, "r") as f:
        pair = f.readline().rstrip("\n").split(" ")

    return OrderedPair(float(pair[0]), float(pair[1]))

def read_center_of_gravity_tests(directory):
    input_files = sorted(os.listdir(directory + "/Input"))
    output_files = sorted(os.listdir(directory + "/Output"))

    # ensure same number of input and output files
    if len(output_files) != len(input_files):
        raise ValueError("Error: number of input and output files do not match!")

    tests = []
    for input_file, output_file in zip(input_files, output_files):
        stars = read_center_of_gravity_sample(directory + "/Input/" + input_file)
        result = read_ordered_pair_from_file(directory + "/Output/" + output_file)
        tests.append(Center_Of_Gravity_Test(stars=stars, result=result))
    return tests

def test_center_of_gravity():
    tests = read_center_of_gravity_tests("center_of_gravity_tests")
    for test in tests:
        our_answer = center_of_gravity(*test.stars)
        # round answer since test result floats are given to 2 decimals
        our_answer.x = round(our_answer.x, 2)
        our_answer.y = round(our_answer.y, 2)

        assert our_answer == test.result, (f"center_of_gravity({test.stars}) = {our_answer}, want {test.result}")


def read_create_children_sample(file):
    with open(file, "r") as f:
        quad = Quadrant()
        quad_attributes = f.readline().rstrip("\n").split(" ")
        quad.width = float(quad_attributes[0])
        quad.x = float(quad_attributes[1])
        quad.y = float(quad_attributes[2])
    return Node(quad)

def read_create_children_result(file):
    children = []
    with open(file, "r") as f:
        for line in f:
            # The each line is child quadrant's attributes
            quad = Quadrant()
            quad_attributes = line.rstrip("\n").split(" ")
            quad.width = float(quad_attributes[0])
            quad.x = float(quad_attributes[1])
            quad.y = float(quad_attributes[2])
            children.append(Node(quad))
    return children

def read_create_children_tests(directory):
    input_files = sorted(os.listdir(directory + "/Input"))
    output_files = sorted(os.listdir(directory + "/Output"))

    # ensure same number of input and output files
    if len(output_files) != len(input_files):
        raise ValueError("Error: number of input and output files do not match!")

    tests = []
    for input_file, output_file in zip(input_files, output_files):
        node = read_create_children_sample(directory + "/Input/" + input_file)
        result = read_create_children_result(directory + "/Output/" + output_file)
        tests.append(Create_Children_Test(node=node, result=result))
    return tests

def test_create_children():
    tests = read_create_children_tests("create_children_tests")
    for test in tests:
        test.node.create_children()

        assert test.node.children == test.result, (f"create_children({test.node}) = {test.node.children}, want {test.result}")


def read_find_child_sample(file):
    with open(file, "r") as f:
        # first line is star to find child for
        star = Star()
        star_attributes = f.readline().rstrip("\n").split(" ")
        star.position = OrderedPair(float(star_attributes[0]), float(star_attributes[1]))
        star.mass = float(star_attributes[2])

        # 2nd line is root node, use create children to generate children
        quad = Quadrant()
        quad_attributes = f.readline().rstrip("\n").split(" ")
        quad.width = float(quad_attributes[0])
        quad.x = float(quad_attributes[1])
        quad.y = float(quad_attributes[2])
        
        root = Node(quad)
        root.create_children()
    return star, root

def read_find_child_result(file):
    with open(file, "r") as f:
        quad = Quadrant()
        quad_attributes = f.readline().rstrip("\n").split(" ")
        quad.width = float(quad_attributes[0])
        quad.x = float(quad_attributes[1])
        quad.y = float(quad_attributes[2])
    return Node(quad)

def read_find_child_tests(directory):
    input_files = sorted(os.listdir(directory + "/Input"))
    output_files = sorted(os.listdir(directory + "/Output"))

    # ensure same number of input and output files
    if len(output_files) != len(input_files):
        raise ValueError("Error: number of input and output files do not match!")

    tests = []
    for input_file, output_file in zip(input_files, output_files):
        star, node = read_find_child_sample(directory + "/Input/" + input_file)
        result = read_find_child_result(directory + "/Output/" + output_file)
        tests.append(Find_Child_Test(star=star, node=node, result=result))
    return tests

def test_find_child():
    tests = read_find_child_tests("find_child_tests")
    for test in tests:
        node = test.node.find_child(test.star)

        assert node == test.result, (f"find_child({test.node, test.star}) = {node}, want {test.result}")


def read_node_insert_sample(file):
    nodes = {}
    root = None
    inserted_star = None

    with open(file, "r") as f:
        # read in each node in quadtree
        for line in f:
            line = line.strip()
            # using empty line to seperate nodes from tree structure and inserted star
            if not line:
                break

            quad = Quadrant()
            star = None
            node_attributes = line.split(":")
            label = node_attributes[0]

            sector_star = node_attributes[1].strip().split("|")
            sector_attributes = sector_star[0].strip().split(" ")
            quad.x = float(sector_attributes[0])
            quad.y = float(sector_attributes[1])
            quad.width = float(sector_attributes[2])
            
            star_attributes = sector_star[1].strip().split(" ")
            if star_attributes[0] != "None":
                star = Star()
                star.position = OrderedPair(float(star_attributes[0]), float(star_attributes[1]))
                star.mass = float(star_attributes[2])

            # use label to add node to dict for later
            nodes[label] = Node(sector=quad, star=star)

        # star to insert
        inserted_star_attributes = f.readline().strip().split(" ")
        inserted_star = Star()
        inserted_star.position = OrderedPair(float(inserted_star_attributes[0]), 
                                                  float(inserted_star_attributes[1]))
        inserted_star.mass = float(inserted_star_attributes[2])

        # read in node relationships, first node is root
        for line in f:
            tree_attributes = line.strip().split("|")
            parent = nodes.get(tree_attributes[0].strip())
            children_labels = tree_attributes[1].strip().split(" ")

            # if root isn't set, set root
            if root == None:
                root = parent

            parent.children = []
            # use dict to append child nodes to parent children list
            for child_label in children_labels:
                parent.children.append(nodes.get(child_label))

        # set root if there's only a single node
        if len(nodes) == 1:
            root = list(nodes.values())[0]

    return root, inserted_star

def read_node_insert_result(file):
    nodes = {}
    root = None

    with open(file, "r") as f:
        # read in each node in quadtree
        for line in f:
            line = line.strip()
            # using empty line to seperate nodes from tree structure
            if not line:
                break

            quad = Quadrant()
            star = None
            node_attributes = line.split(":")
            label = node_attributes[0]

            sector_star = node_attributes[1].strip().split("|")
            sector_attributes = sector_star[0].strip().split(" ")
            quad.x = float(sector_attributes[0])
            quad.y = float(sector_attributes[1])
            quad.width = float(sector_attributes[2])
            
            star_attributes = sector_star[1].strip().split(" ")
            if star_attributes[0] != "None":
                star = Star()
                star.position = OrderedPair(float(star_attributes[0]), float(star_attributes[1]))
                star.mass = float(star_attributes[2])

            # use label to add node to dict for later
            nodes[label] = Node(sector=quad, star=star)

        # read in node relationships, first node is root
        for line in f:
            tree_attributes = line.strip().split("|")
            parent = nodes.get(tree_attributes[0].strip())
            children_labels = tree_attributes[1].strip().split(" ")

            # if root isn't set, set root
            if root == None:
                root = parent
            parent.children = []
            # use dict to append child nodes to parent children list
            for child_label in children_labels:
                parent.children.append(nodes.get(child_label))     

        # set root if there's only a single node
        if len(nodes) == 1:
            root = list(nodes.values())[0]  

    return root

def read_node_insert_tests(directory):
    input_files = sorted(os.listdir(directory + "/Input"))
    output_files = sorted(os.listdir(directory + "/Output"))

    # ensure same number of input and output files
    if len(output_files) != len(input_files):
        raise ValueError("Error: number of input and output files do not match!")

    tests = []
    for input_file, output_file in zip(input_files, output_files):
        root, star = read_node_insert_sample(directory + "/Input/" + input_file)
        result = read_node_insert_result(directory + "/Output/" + output_file)
        tests.append(Node_Insert_Test(root=root, star=star, result=result))
    return tests

def test_node_insert():
    tests = read_node_insert_tests("node_insert_tests")
    
    for test in tests:
        inital_root = copy.deepcopy(test.root)
        test.root.insert(test.star)
        assert test.root == test.result, (f"node_insert({test.star, inital_root}) = {test.root}, want {test.result}")


def read_generate_quadtree_sample(file):
    stars = []
    with open(file, "r") as f:
        width = float(f.readline().strip())
        for line in f:
            # The each line is a star's attributes
            star = Star()
            star_attributes = line.rstrip("\n").split(" ")
            star.position = OrderedPair(float(star_attributes[0]), float(star_attributes[1]))
            star.mass = float(star_attributes[2])
            stars.append(star)
    return Universe(width=width, stars=stars)

def read_generate_quadtree_result(file):
    node_stars = {}
    root = None

    with open(file, "r") as f:
        width = float(f.readline().strip())
        # read in each node in quadtree
        for line in f:
            line = line.strip()
            # using empty line to seperate nodes from tree structure
            if not line:
                break

            star = None
            label_star = line.split(":")
            label = label_star[0]
            
            star_attributes = label_star[1].strip().split(" ")
            if star_attributes[0] != "None":
                star = Star()
                star.position = OrderedPair(float(star_attributes[0]), float(star_attributes[1]))
                star.mass = float(star_attributes[2])

            # use label to add node to dict for later
            node_stars[label] = star

        nodes = {}
        # read in node relationships, first node is root
        for line in f:
            tree_attributes = line.strip().split("|")
            parent_label = tree_attributes[0].strip()
            parent_star = node_stars.get(parent_label)
            children_labels = tree_attributes[1].strip().split(" ")

            # if root isn't set, set root
            if root == None:
                quad = Quadrant(width=width)
                root = Node(sector=quad, star=parent_star)
                nodes[parent_label] = root
            parent = nodes[parent_label]
            parent.create_children()

            # use dict to append child nodes to parent children list
            for i in range(len(parent.children)):
                parent.children[i].star = node_stars.get(children_labels[i])
                nodes[children_labels[i]] = parent.children[i]

        # set root if there's only a single node
        if len(nodes) == 1:
            root = list(nodes.values())[0]  
    
        return QuadTree(root=root)

def read_generate_quadtree_tests(directory):
    input_files = sorted(os.listdir(directory + "/Input"))
    output_files = sorted(os.listdir(directory + "/Output"))

    # ensure same number of input and output files
    if len(output_files) != len(input_files):
        raise ValueError("Error: number of input and output files do not match!")

    tests = []
    for input_file, output_file in zip(input_files, output_files):
        universe = read_generate_quadtree_sample(directory + "/Input/" + input_file)
        result = read_generate_quadtree_result(directory + "/Output/" + output_file)
        tests.append(Generate_QuadTree_Test(universe=universe, result=result))
    return tests

def test_generate_quadtree():
    tests = read_generate_quadtree_tests("generate_quadtree_tests")
    for test in tests:
        tree = generate_quadtree(test.universe)

        assert tree == test.result, (f"generate_quadtree({test.universe}) = {tree}, want {test.result}")


def main():
    test_center_of_gravity()
    test_create_children()
    test_find_child()
    test_node_insert()
    test_generate_quadtree()


if __name__ == "__main__":
    main()
