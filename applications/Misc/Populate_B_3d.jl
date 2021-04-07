module Populate_B_3d

using MATLAB

function init(maxlevel::Integer,folder_with_elements::String,folder::String)

MatlabPopulate() = eval_string(string("Populate_3d_MATLAB(", maxlevel, ")"))
println("###### POPULATING ########")
MatlabPopulate()

end
















end
