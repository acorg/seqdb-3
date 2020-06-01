
function main()
{
    console.log("compare_sequences_data", compare_sequences_data);
    show_full_sequences();
}

// --------------------------------------------------------------------------------

function show_full_sequences()
{
    const add_seqence = function(tr, seq_s) {
        const seq = [...seq_s];
        seq.forEach(function(aa, pos0) {
            const aa_td = document.createElement("td");
            aa_td.classList.add(`aa${aa}`);
            aa_td.classList.add("aa");
            if (pos0 % 10 == 0)
                aa_td.classList.add("sep-left-one");
            else if (pos0 % 10 == 5)
                aa_td.classList.add("sep-left-six");
            aa_td.innerHTML = aa;
            tr.appendChild(aa_td);
        });
    };

    const add_ruler = function() {
        const tr = document.createElement("tr");
        tr.classList.add("aa-ruler");
        tr.innerHTML = `<td colspan="2"></td>`;
        for (let pos1 = 1; pos1 < 550; ++pos1) {
            const td = document.createElement("td");
            if (pos1 % 10 == 1) {
                td.classList.add("sep-left-one");
                let pos_s = "" + pos1;
                td.innerHTML = pos_s;
                td.setAttribute("colspan", pos_s.length);
                pos1 += pos_s.length - 1;
            }
            else if (pos1 % 10 == 6) {
                td.classList.add("sep-left-six");
                td.innerHTML = ".";
            }
            else
                td.innerHTML = ".";
            tr.appendChild(td);
        }
        return tr;
    };

    const tab1 = document.createElement("table");
    for (let group of compare_sequences_data.groups) {
        tab1.appendChild(add_ruler());
        group.seq.forEach(function(id_seq, index) {
            const tr = document.createElement("tr");
            if (index == 0)
                tr.innerHTML = `<td class="group-name" rowspan="${group.seq.length}">${group.name}</td><td class="seq-id">${id_seq.id}</td>`;
            else
                tr.innerHTML = `<td class="seq-id">${id_seq.id}</td>`;
            add_seqence(tr, id_seq.seq);
            tab1.appendChild(tr);
        });
    }
    tab1.appendChild(add_ruler());
    document.querySelector("#full-sequences").appendChild(tab1);
}

// --------------------------------------------------------------------------------

// --------------------------------------------------------------------------------

document.addEventListener("DOMContentLoaded", main);
