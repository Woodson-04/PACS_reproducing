param(
  [string]$InputMarkdown = "PACS_GSE157079_UMAP_progress_report_FINAL.md",
  [string]$OutputHtml = "PACS_GSE157079_UMAP_progress_report_FINAL.html"
)

function HtmlEncode([string]$s) {
  if ($null -eq $s) { return "" }
  return [System.Net.WebUtility]::HtmlEncode($s)
}

function InlineFormat([string]$s) {
  $encoded = HtmlEncode $s
  $encoded = [regex]::Replace($encoded, '!\[([^\]]*)\]\(([^)]+)\)', '<img alt="$1" src="$2" />')
  $encoded = [regex]::Replace($encoded, '\[([^\]]+)\]\(([^)]+)\)', '<a href="$2">$1</a>')
  $encoded = [regex]::Replace($encoded, '\*\*([^*]+)\*\*', '<strong>$1</strong>')
  $encoded = [regex]::Replace($encoded, '`([^`]+)`', '<code>$1</code>')
  return $encoded
}

function IsSeparatorRow([string]$line) {
  return $line -match '^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)+\|?\s*$'
}

function SplitTableRow([string]$line) {
  $trimmed = $line.Trim()
  if ($trimmed.StartsWith("|")) { $trimmed = $trimmed.Substring(1) }
  if ($trimmed.EndsWith("|")) { $trimmed = $trimmed.Substring(0, $trimmed.Length - 1) }
  return $trimmed -split '\s*\|\s*'
}

$lines = Get-Content -Path $InputMarkdown -Encoding UTF8
$out = New-Object System.Collections.Generic.List[string]
$out.Add('<!doctype html>')
$out.Add('<html lang="zh-CN">')
$out.Add('<head>')
$out.Add('<meta charset="utf-8" />')
$out.Add('<meta name="viewport" content="width=device-width, initial-scale=1" />')
$out.Add('<title>PACS 复现项目阶段报告</title>')
$out.Add('<link rel="stylesheet" href="report_style.css" />')
$out.Add('<script>')
$out.Add('window.MathJax = { tex: { inlineMath: [["\\(","\\)"]], displayMath: [["$$","$$"]] }, svg: { fontCache: "global" } };')
$out.Add('</script>')
$out.Add('<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>')
$out.Add('</head>')
$out.Add('<body>')

$inCode = $false
$codeBuffer = New-Object System.Collections.Generic.List[string]
$inUl = $false
$inOl = $false
$inTable = $false

for ($i = 0; $i -lt $lines.Count; $i++) {
  $line = $lines[$i]

  if ($line -match '^```') {
    if (-not $inCode) {
      $inCode = $true
      $codeBuffer.Clear()
    } else {
      $inCode = $false
      $out.Add('<pre><code>' + (HtmlEncode (($codeBuffer -join "`n"))) + '</code></pre>')
    }
    continue
  }

  if ($inCode) {
    $codeBuffer.Add($line)
    continue
  }

  if ($inTable -and ($line -notmatch '^\s*\|')) {
    $out.Add('</tbody></table>')
    $inTable = $false
  }
  if ($inUl -and ($line -notmatch '^\s*-\s+')) {
    $out.Add('</ul>')
    $inUl = $false
  }
  if ($inOl -and ($line -notmatch '^\s*\d+\.\s+')) {
    $out.Add('</ol>')
    $inOl = $false
  }

  if ([string]::IsNullOrWhiteSpace($line)) {
    continue
  }

  if ($line -match '^---+\s*$') {
    $out.Add('<hr />')
    continue
  }

  if ($line -match '^(#{1,6})\s+(.+)$') {
    $level = $matches[1].Length
    $text = InlineFormat $matches[2]
    $out.Add("<h$level>$text</h$level>")
    continue
  }

  if ($line -match '^\s*>\s*(.+)$') {
    $out.Add('<blockquote>' + (InlineFormat $matches[1]) + '</blockquote>')
    continue
  }

  if ($line -match '^\s*\|') {
    if (($i + 1) -lt $lines.Count -and (IsSeparatorRow $lines[$i + 1])) {
      $headers = SplitTableRow $line
      $out.Add('<table><thead><tr>')
      foreach ($h in $headers) { $out.Add('<th>' + (InlineFormat $h.Trim()) + '</th>') }
      $out.Add('</tr></thead><tbody>')
      $inTable = $true
      $i++
      continue
    } elseif ($inTable) {
      $cells = SplitTableRow $line
      $out.Add('<tr>')
      foreach ($c in $cells) { $out.Add('<td>' + (InlineFormat $c.Trim()) + '</td>') }
      $out.Add('</tr>')
      continue
    }
  }

  if ($line -match '^\s*-\s+(.+)$') {
    if (-not $inUl) {
      $out.Add('<ul>')
      $inUl = $true
    }
    $out.Add('<li>' + (InlineFormat $matches[1]) + '</li>')
    continue
  }

  if ($line -match '^\s*(\d+)\.\s+(.+)$') {
    if (-not $inOl) {
      $out.Add('<ol>')
      $inOl = $true
    }
    $out.Add('<li>' + (InlineFormat $matches[2]) + '</li>')
    continue
  }

  if ($line -match '^\s*!\[([^\]]*)\]\(([^)]+)\)\s*$') {
    $out.Add('<figure><img alt="' + (HtmlEncode $matches[1]) + '" src="' + (HtmlEncode $matches[2]) + '" /></figure>')
    continue
  }

  $out.Add('<p>' + (InlineFormat $line) + '</p>')
}

if ($inTable) { $out.Add('</tbody></table>') }
if ($inUl) { $out.Add('</ul>') }
if ($inOl) { $out.Add('</ol>') }
if ($inCode) { $out.Add('<pre><code>' + (HtmlEncode (($codeBuffer -join "`n"))) + '</code></pre>') }

$out.Add('</body>')
$out.Add('</html>')

$htmlText = $out -join "`n"
[System.IO.File]::WriteAllText((Join-Path (Get-Location) $OutputHtml), $htmlText, [System.Text.Encoding]::UTF8)
