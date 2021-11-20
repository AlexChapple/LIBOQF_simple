"""

Modifies the Little_Boxes_Simple_3D.f08 file and runs a shell script each time to run the 
code in the terminal 

"""

import os 

# makes the list of increments we want to make 
end_time_list = range(2,102,2)

# Iterate over the different end_times and generate the data 
for index in range(len(end_time_list)):

    # Edit the Little_Boxes_Simple_3D.f08 file 

    file = open("Little_Boxes_Simple_3D.f08", "r")
    list_of_lines = file.readlines()
    list_of_lines[14] = "    integer, parameter :: end_time = " + str(end_time_list[index]) + "d0\n"

    file = open("Little_Boxes_Simple_3D.f08", "w")
    file.writelines(list_of_lines)
    file.close()

    # Run the Little_Boxes_Simple_3D.f08 in the terminal 
    os.system("gfortran-9 Little_Boxes_Simple_3D.f08 -o LIBOQF_3D.out")
    os.system("./LIBOQF_3D.out")

    print(str(index + 1) + "/" + str(len(end_time_list)) + " end_time simulations completed")




