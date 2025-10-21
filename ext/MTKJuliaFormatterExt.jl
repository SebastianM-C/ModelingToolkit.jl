module MTKJuliaFormatterExt

using ModelingToolkit
using JuliaFormatter

"""
    readable_code(expr)

Format an expression for readable output using JuliaFormatter.

This function is only available when JuliaFormatter is loaded. Without JuliaFormatter,
the expression will be returned as a string without additional formatting.
"""
function ModelingToolkit.readable_code(expr)
    expr = Base.remove_linenums!(ModelingToolkit._readable_code(expr))
    ModelingToolkit.rec_remove_macro_linenums!(expr)
    JuliaFormatter.format_text(string(expr), JuliaFormatter.SciMLStyle())
end

end
