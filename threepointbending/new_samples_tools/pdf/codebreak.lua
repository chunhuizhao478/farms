-- codebreak.lua: emit inline Code as \texttt{...} with break opportunities so
-- long repo paths, slugs and hashes wrap inside table cells and prose instead
-- of running off the page.  LaTeX output only.
--
--   * after '/', '_', '-', '.', '=', ':' and ','  ->  \allowbreak{}   (free break)
--   * after any other character of a token longer than 16 chars
--                                               ->  \penalty5000{}   (last resort)
--   * ordinary spaces stay breakable glue (pandoc's own \texttt uses control
--     spaces, which never break)
--   * the pipe-table escape \| becomes a plain |
local special = {
  ['{'] = '\\{', ['}'] = '\\}', ['#'] = '\\#', ['%'] = '\\%', ['&'] = '\\&',
  ['$'] = '\\$', ['^'] = '\\^{}', ['~'] = '\\textasciitilde{}',
  ['\\'] = '\\textbackslash{}', ['_'] = '\\_',
}
local soft = { ['/'] = true, ['_'] = true, ['-'] = true, ['.'] = true,
               ['='] = true, [':'] = true, [','] = true }

local function esc(s)
  s = s:gsub('\\|', '|')
  local n = #s
  local long = n > 16
  local out = {}
  for i = 1, n do
    local c = s:sub(i, i)
    -- keep command-line option dashes glued: no break inside '--' or after a dash that
    -- starts an option ('-m', '--blend-km'); ordinary in-word dashes stay breakable.
    -- TeX may break after ANY literal hyphen, so an option dash is boxed.
    local nxt, prv = s:sub(i + 1, i + 1), (i > 1) and s:sub(i - 1, i - 1) or ' '
    local optdash = (c == '-') and (nxt == '-' or prv == ' ' or prv == '-')
    out[#out + 1] = optdash and '\\mbox{-}' or (special[c] or c)
    if i < n and c ~= ' ' and not optdash then
      if soft[c] then
        out[#out + 1] = '\\allowbreak{}'
      elseif long then
        out[#out + 1] = '\\penalty5000{}'
      end
    end
  end
  return table.concat(out)
end

function Code(el)
  if FORMAT:match('latex') then
    return pandoc.RawInline('latex', '\\texttt{' .. esc(el.text) .. '}')
  end
end
