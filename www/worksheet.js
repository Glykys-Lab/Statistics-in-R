// Browser bridge for worksheet edits, navigator clicks and saved analysis defaults.
// Send events to Shiny; validation and calculations belong to the R server.
// A worksheet revision identifies the rendered table, and a nonce makes repeated
// actions observable even when their other payload fields have not changed.
function saveWorksheetCell(event) {
  const cell = event.target.closest('.worksheet-cell');
  if (!cell || !window.Shiny) return;
  Shiny.setInputValue('cell_edit', {
    row: Number(cell.dataset.row), column: Number(cell.dataset.column),
    value: cell.value, revision: Number(cell.dataset.revision), nonce: Date.now()
  }, {priority: 'event'});
}
// Capture before Shiny's ordinary text-input bindings handle the same event.
document.addEventListener('change', saveWorksheetCell, true);
document.addEventListener('blur', saveWorksheetCell, true);
document.addEventListener('click', function(event) {
  // Delegate clicks because Shiny replaces navigator buttons as the project changes.
  const link = event.target.closest('.project-link');
  if (link && window.Shiny) Shiny.setInputValue('project_open', {
    kind: link.dataset.kind, key: link.dataset.key, nonce: Date.now()
  }, {priority: 'event'});
}, true);
// Keyboard movement within the worksheet keeps the familiar spreadsheet order.
document.addEventListener('keydown', function (event) {
  const cell = event.target.closest('.worksheet-cell');
  if (!cell || !['ArrowUp', 'ArrowDown', 'Enter'].includes(event.key)) return;
  event.preventDefault();
  const row = Number(cell.dataset.row) + (event.key === 'ArrowUp' ? -1 : 1);
  const next = document.querySelector('.worksheet-cell[data-row="' + row + '"][data-column="' + cell.dataset.column + '"]');
  if (next) { next.focus(); next.select(); }
});

// Browser-local defaults contain analysis choices, never worksheet measurements.
function initializeAnalysisDefaults() {
  // Bind only after Shiny is available and restore storage after its connection
  // is ready. The server validates these saved choices before using them.
  if (!window.Shiny || !window.jQuery) {
    window.setTimeout(initializeAnalysisDefaults, 50);
    return;
  }
  Shiny.addCustomMessageHandler('analysis-defaults', function (message) {
    try {
      if (message.clear) localStorage.removeItem('labstats.analysis.defaults.v1');
      else localStorage.setItem('labstats.analysis.defaults.v1', JSON.stringify(message.values));
    } catch (_) { Shiny.setInputValue('defaults_storage_error', Date.now(), {priority: 'event'}); }
  });
  function restoreAnalysisDefaults() {
    try {
      const values = JSON.parse(localStorage.getItem('labstats.analysis.defaults.v1'));
      if (values && typeof values === 'object' && !Array.isArray(values)) {
        Shiny.setInputValue('restore_analysis_defaults', values, {priority: 'event'});
      }
    } catch (_) {
      try { localStorage.removeItem('labstats.analysis.defaults.v1'); } catch (_) {}
    }
  }
  if (Shiny.shinyapp && Shiny.shinyapp.isConnected()) restoreAnalysisDefaults();
  else jQuery(document).one('shiny:connected', restoreAnalysisDefaults);
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initializeAnalysisDefaults);
else initializeAnalysisDefaults();
