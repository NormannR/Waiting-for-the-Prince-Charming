module Macros

export @def, @load_struct, @load_struct_cat, @save_struct, @save_struct_index, @load_index, @load_index_cat, init!, @load_vec_into_struct

"""
	@def name definition

Defines the macro `name`, which writes the code in `definition`
"""
macro def(name, definition)
  return quote
      macro $(esc(name))()
          esc($(Expr(:quote, definition)))
      end
  end
end

"""
	@load_struct_cat name s v1 v2 ...

Load variables `v1`, `v2`, ... from structure `s` and rename them as suffix `v1_name`, `v2_name`, ...
"""
macro load_struct_cat(name, d, vars...)
	expr = Expr(:block)
	for v in vars
		push!(expr.args, :($(Symbol("$(v)_$(name)")) = $d.$v))
	end
	return esc(expr)
end

"""
	@load_index d i v1 v2 ...

Loads variables `v1`, `v2`, ... from structure `d` from the index `i` of the associated fields.
"""
macro load_index(d, i, vars...)
	expr = Expr(:block)
	for v in vars
		push!(expr.args, :($v = ($d.$v)[$i])) 
	end
	return esc(expr)
end

"""
	@load_index_cat d i cat v1 v2 ...

Loads variables `v1_cat`, `v2_cat`, ... from structure `d` from the index `i` of the associated fields `v1`, `v2`, ... .
"""
macro load_index_cat(d, i, cat, vars...)
	expr = Expr(:block)
	for v in vars
		push!(expr.args, :($(Symbol("$(v)_$(cat)")) = ($d.$v)[$i])) 
	end
	return esc(expr)
end

"""
	@load_struct d v1 ...

Loads variables `v1`, ... from the structure `d`
"""
macro load_struct(d, vars...)
	expr = Expr(:block)
	for v in vars
		push!(expr.args, :($v = $d.$v))
	end
	return esc(expr)
end

"""
	@save_struct d v1 ...

Saves local variables `v1`, ... in the associated slot `d.v1` of structure `d`
"""
macro save_struct(d, vars...)
	expr = Expr(:block)
	for v in vars
		push!(expr.args, :($d.$v = $v))
	end
	return esc(expr)
end

"""
	@save_struct_index d i v1 ...

Saves local variables `v1`, ... at the index `i` of the associated slot `(d.v1)[i]` of structure `d`
"""
macro save_struct_index(d, i, vars...)
	expr = Expr(:block)
	for v in vars
		push!(expr.args, :(($d.$v)[$i] = $v))
	end
	return esc(expr)
end

"""
	@init_struct s v

Initializes the properties of structure `s` to value `v`
"""
macro init_struct(s, v)
	expr = Expr(:block)
	for v in vars
		push!(expr.args, :($s.$v = $v))
	end
	return expr
end

"""
	init!(s::S)
"""
@generated function init!(s::S, N) where S
	expr = Expr(:block)
	for v in fieldnames(S)
		push!(expr.args, :(s.$v = zeros(N)))
	end
	return expr
end

"""
	@load_vec_into_struct vec s v1 ...

Assign values `vec[1]`, `vec[2]`, ... to  fields `v1`, `v2`, ... of structure `s`
"""
macro load_vec_into_struct(vec, s, vars...)
	expr = Expr(:block)
	for (k,v) in enumerate(vars)
		push!(expr.args, :($s.$v = $vec[$k]))
	end
	return esc(expr)
end

end
