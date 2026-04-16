local function pick_delim(s)
  local delimiters = {"|", "/", "!", "+", ";", "~", "^"}
  for _, d in ipairs(delimiters) do
    if not s:find(d, 1, true) then
      return d
    end
  end
  return "|"
end

local function escape_latex_texttt(s)
  s = s:gsub("\\", "\\textbackslash{}")
  s = s:gsub("([{}%%#$&_])", "\\%1")
  s = s:gsub("%^", "\\textasciicircum{}")
  s = s:gsub("~", "\\textasciitilde{}")
  return s
end

local function heading_code(code)
  local s = code.text
  local tt = escape_latex_texttt(s)
  local plain = s:gsub("[{}\\]", "")
  return pandoc.RawInline(
    "latex",
    "\\texorpdfstring{\\texttt{" .. tt .. "}}{" .. plain .. "}"
  )
end

local function body_code(code)
  local s = code.text
  local delim = pick_delim(s)
  return pandoc.RawInline(
    "latex",
    "\\Verb[breaklines=true,breakanywhere=true,breakanywhereinlinestretch=1em,bgcolor=codebg]" ..
      delim .. s .. delim
  )
end

function Pandoc(doc)
  for i, blk in ipairs(doc.blocks) do
    if blk.t == "Header" then
      doc.blocks[i] = pandoc.walk_block(blk, {
        Code = heading_code
      })
    else
      doc.blocks[i] = pandoc.walk_block(blk, {
        Code = body_code
      })
    end
  end
  return doc
end